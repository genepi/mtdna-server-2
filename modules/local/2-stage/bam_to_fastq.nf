process BAM_TO_FASTQ {

  input:
    path bamFile

  output:
    tuple val("${bamFile.baseName}"), path("*fastq"), emit: fastq_ch

	"""
    bedtools bamtofastq -i $bamFile -fq ${bamFile.baseName}.fastq
	"""
}