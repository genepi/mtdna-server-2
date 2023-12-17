process INDEX {
	input:
	path reference

	output:
	path "*.{dict,fai}", emit: fasta_index_ch

	"""
    samtools faidx $reference
    samtools dict $reference -o ${reference.baseName}.dict
	"""
}