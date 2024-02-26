process FILTER_VARIANTS {

    //publishDir "${params.output}", mode: 'copy'

    input:
    tuple path(vcf_file), val(method)

    output:
    path("${vcf_file.baseName}.${method}.filtered.txt"), emit: combined_methods_ch

    script:
    def vcf_name = "${vcf_file}".replaceAll('.vcf.gz', '')

    """
    echo -e "ID\tFilter\tPos\tRef\tVariant\tVariantLevel\tCoverage\tType" > ${vcf_file.baseName}.${method}.txt
    
    if [[ ${method} == "mutserve" ]]
    then
        bcftools query -f '${vcf_name}.bam\t%FILTER\t%POS\t%REF\t%ALT\t[%AF\t%DP\t%GT]\n' ${vcf_file} >> ${vcf_file.baseName}.${method}.txt
        awk -F'\t' 'NR == 1 || (length(\$4) == 1 && length(\$5) == 1)' ${vcf_file.baseName}.mutserve.txt > ${vcf_file.baseName}.${method}.filtered.txt
    else
        bcftools query -f '${vcf_name}.bam\t%FILTER\t%POS\t%REF\t%ALT\t[%AF\t%DP]\tINDEL\n'  ${vcf_file} >> ${vcf_file.baseName}.${method}.txt
        awk -F'\t' 'NR == 1 || ((length(\$4) > 1 || length(\$5) > 1) && length(\$4) != length(\$5))' ${vcf_file.baseName}.${method}.txt > ${vcf_file.baseName}.${method}.filtered.txt
    fi
    """
}

