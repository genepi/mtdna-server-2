process QUALITY_CONTROL {
    
    publishDir "${params.output_reports}/multiqc", mode: "copy", pattern: '*.html'
    
    input:
    path excluded_samples
    path zip
    
    output:
	path "*.html"

	"""
	multiqc . --filename index.html
	"""
}