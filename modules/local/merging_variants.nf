process MERGING_VARIANTS {

    input:
    path variants_txt
    val mode

    output:
    path("variants.txt"), emit: txt_summarized_ch

    """
    csvtk concat \
        -t ${variants_txt} \
        -T -o variants.concat.txt \
        --num-cpus ${task.cpus}
    
    csvtk sort \
        -t variants.concat.txt \
        -k ID:N -k Pos:n -k Ref:N -k Variant:N -k Type:nr \
        -T -o variants.sorted.txt \
        --num-cpus ${task.cpus}

    if [[ ${mode} == "fusion" ]]
    then
        java -jar /opt/VariantMerger.jar \
            variants.sorted.txt \
            --output variants.txt
    else
        mv variants.sorted.txt variants.txt
    fi
    """
}