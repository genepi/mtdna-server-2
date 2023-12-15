process QUALITY_CONTROL {
    publishDir "${params.output}", mode: "copy", pattern: '*.html'
    input:
    path excluded_samples
    path zip
    output:
	path "*.html"

	"""
	 multiqc . 
	"""
}