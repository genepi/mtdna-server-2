process COVERAGE_ESTIMATION {

    input:
    path variants_txt

    output:
    path("variants.txt"), emit: variants_verified_ch

    """
    java -jar /opt/CoverageEstimation.jar \
        --input ${variants_txt} \
        --output variants.verified.txt
    mv variants.verified.txt variants.txt     
    """
}