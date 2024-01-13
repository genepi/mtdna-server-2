process PAIRWISE_ALIGNMENT {
	input:
	path fasta1
	path fasta2

	output:
	path "*.fasta", emit: aligned

	"""
    water \
	-auto \
	-stdout \
	-asequence /home/seb/Desktop/rCRS.fasta \
	-bsequence /home/seb/Desktop/chrM.30.fa \
	-datafile EDNAFULL \
	-gapopen 10.0 \
	-gapextend 0.5 \
	-aformat3 fasta \
	-snucleotide1 \
	-snucleotide2
	"""
}