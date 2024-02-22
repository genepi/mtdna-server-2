process STAGE_2 {

  input:
    tuple val(baseName), path(fastq)
    path ref_fasta

  output:
    path "*.2stage.bam", emit: realigned_ch

	"""
  bwa index "${ref_fasta}"
	bwa mem ${ref_fasta} ${fastq} | samtools sort -o ${baseName}.2stage.bam -
	"""
}