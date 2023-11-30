process REPORT {

    publishDir "${params.output}/reports", mode: 'copy'

    input:
    path variants_ch

    output:
    path ("report.html"), emit: report_ch

    """
    java -jar /opt/mutserve/mutserve.jar report --input ${variants_ch} --output report.html
    """
}