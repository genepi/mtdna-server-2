process CONTAMINATION_DETECTION {

    input:
    path variants_ch

    output:
    path ("haplocheck.html"), emit: contamination_report_ch
    """
    java -jar /opt/haplocheck/haplocheck.jar --out haplocheck.txt --raw ${variants_ch}
    """
}