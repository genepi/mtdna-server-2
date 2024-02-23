process SUMMARIZE_VARIANTS {

    publishDir "${params.output}", mode: 'copy'

    input:
    path variants_txt

    output:
    path("variants.fixed.txt"), emit: txt_summarized_ch

    """
    csvtk concat -t ${variants_txt} -T -o variants.concat.txt
    
    csvtk sort -t variants.concat.txt -k ID:N -k Pos:n -k Type -T -o variants.sorted.txt

    java -jar /opt/VariantMerger.jar \
        variants.sorted.txt \
        --output variants.fixed.txt
    """
}