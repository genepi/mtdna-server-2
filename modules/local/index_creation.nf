process INDEX_CREATION {
	input:
	path reference
	val mtdna_tag

	output:
	path "ref*.{dict,fai}", emit: fasta_index_ch
	path "ref.fasta", emit: ref_ch

	"""
	sed -e "s/^>.*/>$mtdna_tag/" $reference > ref.fasta
    samtools faidx ref.fasta
    	samtools dict ref.fasta \
	    -o ref.dict
	"""
}