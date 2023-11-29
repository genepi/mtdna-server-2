process CONTAMINATION_DETECTION {

    publishDir "${params.output}/reports", mode: 'copy'

    input:
    path variants_ch

    output:
    path ("contamination.html"), emit: contamination_report_ch
    """
    java -jar /opt/haplocheck/haplocheck.jar --out contamination.txt --raw ${variants_ch}
    """
}