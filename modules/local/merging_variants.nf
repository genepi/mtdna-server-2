process MERGING_VARIANTS {

    //publishDir "${params.output}", mode: 'copy'

    input:
    path variants_txt
    val mode
    output:
    path("variants.merged.txt"), emit: txt_summarized_ch

    """
    csvtk concat -t ${variants_txt} -T -o variants.concat.txt
    
    csvtk sort -t variants.concat.txt -k ID:N -k Pos:n -k Type:r -T -o variants.sorted.txt

    if [[ ${mode} == "fusion" ]]
    then
        java -jar /opt/VariantMerger.jar \
            variants.sorted.txt \
            --output variants.merged.txt
    else
        mv variants.sorted.txt variants.merged.txt
    fi
    """
}