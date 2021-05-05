version 1.0

# MobiDL 2.0 - MobiDL 2 is a collection of tools wrapped in WDL to be used in any WDL pipelines.
# Copyright (C) 2021 MoBiDiC
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import "tasks/utilities.wdl" as utilities
import "tasks/GATK4.wdl" as GATK4
import "tasks/sambamba.wdl" as sambamba
import "tasks/star.wdl" as star

workflow GATK_best_practices_RNA {
	meta {
		author: "MoBiDiC"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-04-07"
		source: "https://gatk.broadinstitute.org/hc/en-us/articles/360035531192"
	}

################################################################################
## Inputs

	input {
		File fastqR1
		File? fastqR2

		File refFasta
		File refFai
		File refDict

		File genomeDir

		File? regionOfInterest

		File dbsnp
		File dbsnpIdx

		Array[File] knownSites
		Array[File] knownSitesIdx

		String outputPath
		String subString = "(_S[0-9]+)?(_L[0-9]+)?([_.]?R[12])?(_[0-9]+)?.(fastq|fq)(.gz)?"
		String subStringReplace = ""
		String? name

		Int scatterCount = 4
	}

	String sampleName = if defined(name) then "~{name}" else sub(basename(fastqR1),subString,subStringReplace)

################################################################################

################################################################################
## Alignment

	call star.alignReads as ALIGNREADS {
		input :
			fastqR1 = fastqR1,
			fastqR2 = fastqR2,
			name = name,
			genomeDir = genomeDir,
			outputPath = outputPath + "/Alignment/"
	}

	call sambamba.markdup as MARKDUP {
		input :
			in = ALIGNREADS.bam,
			outputPath = outputPath + "/Alignment/"
	}

	call sambamba.sort as SORT {
		input :
			in = MARKDUP.outputBam,
			outputPath = outputPath + "/Alignment/"
	}

	call GATK4.splitNCigarReads as SPLITNCIGARREADS {
		input :
			in = SORT.outputBam,

			refFasta = refFasta,
			refFai = refFai,
			refDict = refDict,

			outputPath = outputPath + "/Alignment/"
	}

################################################################################

################################################################################
## Base recalibrator & LeftAlign

	call utilities.fai2bed as FAI2BED {
		input :
			in = refFai,
			outputPath = outputPath + "/ROI/"
	}

	call utilities.convertBedToIntervals as BED2INTERVALS {
		input :
			in = select_first([regionOfInterest, FAI2BED.outputFile]),
			outputPath = outputPath + "/ROI/"
	}

	call GATK4.splitIntervals as SPLITINTERVALS {
		input :
			in = BED2INTERVALS.outputFile,

			refFasta = refFasta,
			refFai = refFai,
			refDict = refDict,

			subdivisionMode = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION",

			outputPath = outputPath + "/regionOfInterest/",

			scatterCount = scatterCount
	}

	scatter (intervals in SPLITINTERVALS.splittedIntervals) {
		call GATK4.baseRecalibrator as BASERECALIBRATOR {
			input :
				in = SPLITNCIGARREADS.outputBam,
				bamIdx = SPLITNCIGARREADS.outputIdx,

				intervals = intervals,

				knownSites = knownSites,
				knownSitesIdx = knownSitesIdx,

				refFasta = refFasta,
				refFai = refFai,
				refDict = refDict,

				outputPath = outputPath + "/Alignment/BQSR/"
		}
	}

	call GATK4.gatherBQSRReports as GATHERBQSRREPORTS {
		input :
			in = BASERECALIBRATOR.outputFile,
			outputPath = outputPath + "/Alignment/BQSR/"
	}

	scatter (intervals in SPLITINTERVALS.splittedIntervals) {
		call GATK4.applyBQSR as APPLYBQSR {
			input :
				in = SPLITNCIGARREADS.outputBam,
				bamIdx = SPLITNCIGARREADS.outputIdx,
				bqsrReport = GATHERBQSRREPORTS.outputFile,
				intervals = intervals,

				refFasta = refFasta,
				refFai = refFai,
				refDict = refDict,

				outputPath = outputPath + "/Alignment/BQSR/"
		}

		call GATK4.leftAlignIndels as LEFTALIGNINDELS {
			input :
				in = APPLYBQSR.outputBam,
				bamIdx = APPLYBQSR.outputBai,
				intervals = intervals,

				refFasta = refFasta,
				refFai = refFai,
				refDict = refDict,

				outputPath = outputPath + "/Alignment/BQSR/"
		}
	}

	call GATK4.gatherBamFiles as GATHERBAMFILES {
		input :
			in = LEFTALIGNINDELS.outputBam,
			bamIdx = LEFTALIGNINDELS.outputBai,

			outputPath = outputPath + "/Alignment/"
	}

################################################################################

################################################################################
## VariantCalling

	### HaplotypeCaller
	scatter (intervals in SPLITINTERVALS.splittedIntervals) {
		call GATK4.haplotypeCaller as HAPLOTYPECALLER {
			input :
				in = GATHERBAMFILES.outputBam,
				idx = GATHERBAMFILES.outputBai,

				refFasta = refFasta,
				refFai = refFai,
				refDict = refDict,

				name = sampleName,

				dbsnp = dbsnp,
				dbsnpIdx = dbsnpIdx,

				intervals = intervals,

				outputPath = outputPath + "/Variant-Calling/HaplotypeCaller/"
		}
	}

	call GATK4.gatherVcfFiles as GATHERVCFFILES {
		input :
			in = HAPLOTYPECALLER.outputFile,
			outputPath = outputPath + "/Variant-Calling/HaplotypeCaller/"
	}

	call GATK4.variantFiltration as VARIANTFILTRATION {
		input :
			in = GATHERVCFFILES.outputFile,

			refFasta = refFasta,
			refFai = refFai,
			refDict = refDict,

			clusterSize = 3,
			clusterWindow = 35,

			LowQualByDepth = 2,
			FSStrandBias = 30,

			outputPath = outputPath + "/Variant-Calling/HaplotypeCaller/"
	}

################################################################################

################################################################################
## Outputs

	output {
		File bam = GATHERBAMFILES.outputBam
		File bai = GATHERBAMFILES.outputBai
		File vcfRaw = GATHERVCFFILES.outputFile
		File vcfRawIdx = GATHERVCFFILES.outputFileIdx
		File vcf = VARIANTFILTRATION.outputFile
		File vcfIdx = VARIANTFILTRATION.outputFileIdx
	}

################################################################################

	parameter_meta {
		fastqR1 : {
			description: 'Input file with reads 1 (fastq, fastq.gz, fq, fq.gz).',
			category: 'Required'
		}
		fastqR2 : {
			description: 'Input file with reads 2 (fastq, fastq.gz, fq, fq.gz).',
			category: 'Required'
		}
		refFasta: {
			description: 'Path to the reference file (format: fasta)',
			category: 'Required'
		}
		refFai: {
			description: 'Path to the reference file index (format: fai)',
			category: 'Required'
		}
		refDict: {
			description: 'Path to the reference file dict (format: dict)',
			category: 'Required'
		}
		genomeDir: {
			description: 'Path to the directory where genome files are stored',
			category: 'Required'
		}
		regionOfInterest: {
			description: "Input bed of ROI (default: whole genome)",
			category: "Option"
		}
		knownSites : {
			description: 'Path to the knownSites vcf files (format: vcf.gz)',
			category: 'Required'
		}
		knownSitesIdx : {
			description: 'Path to the knownSites vcf indexes (format: vcf.gz.tbi)',
			category: 'Required'
		}
		dbsnp : {
			description: 'Path to the dbsnp vcf file (format: vcf.gz)',
			category: 'Required'
		}
		dbsnpIdx : {
			description: 'Path to the dbsnp vcf index (format: vcf.gz.tbi)',
			category: 'Required'
		}
		outputPath : {
			description: 'Path where the output repository will be created [default: ./]',
			category: 'Output path/name option'
		}
		subString : {
			description: 'The regexp substring to remove from fastq R1 file to create sampleName, and used by bwa-mem [default: "(_S[0-9]+)?(_L[0-9][0-9][0-9])?(_R[12])?(_[0-9][0-9][0-9])?.(fastq|fq)(.gz)?"]',
			category: 'Output path/name option'
		}
		subStringReplace: {
			description: 'subString replace by this string [default: ""]',
			category: 'Output path/name option'
		}
		name : {
			description: 'The name used as sampleName',
			category: 'Output path/name option'
		}
		scatterCount: {
			description: 'Scatter count: number of output interval files to split into [default: 1]',
			category: 'Parralelisation ption'
		}
	}
}
