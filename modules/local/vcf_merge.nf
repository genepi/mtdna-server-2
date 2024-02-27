process VCF_MERGE {

    input:
    path(variants_ch)
    path(variants_idx_ch)
    val count

    output:
    path ("merged.vcf.gz"), emit: vcf_merged_ch

    """
    if [[ ${count} > 1 ]]
    then
        bcftools merge -Oz -o merged.vcf.gz ${variants_ch} 
    else
        mv ${variants_ch} merged.vcf.gz 
    fi
    """
}