process FILTER_VARIANTS {

    //publishDir "${params.output}", mode: 'copy'

    input:
    tuple path(vcf_file), path(vcf_file_idx), val(method)

    output:
    path("${vcf_file.baseName}.${method}.filtered.txt"), emit: combined_methods_ch

    script:
    def vcf_name = "${vcf_file}".replaceAll('.vcf.gz', '')

    """
    echo -e "ID\tFilter\tPos\tRef\tVariant\tVariantLevel\tCoverage\tGT" \
        > ${vcf_file.baseName}.${method}.txt

    bcftools query \
        -f '${vcf_name}.bam\t%FILTER\t%POS\t%REF\t%ALT\t[%AF\t%DP\t%GT]\n' \
        ${vcf_file} >> ${vcf_file.baseName}.${method}.txt    
    
    if [[ ${method} == "mutserve_fusion" ]]
    then
        awk -F'\t' 'NR == 1 || (length(\$4) == 1 && length(\$5) == 1)' \
            ${vcf_file.baseName}.${method}.txt > ${vcf_file.baseName}.${method}.filtered.tmp.txt

    elif [[ ${method} == "mutect2_fusion" ]]
    then
        awk -F'\t' 'NR == 1 || ((length(\$4) > 1 || length(\$5) > 1) && length(\$4) != length(\$5))' \
            ${vcf_file.baseName}.${method}.txt > ${vcf_file.baseName}.${method}.filtered.tmp.txt
    else 
        mv ${vcf_file.baseName}.${method}.txt ${vcf_file.baseName}.${method}.filtered.tmp.txt  
    fi
    
    ## annotating SNVS and INDELs for reporting
    awk 'BEGIN {OFS="\t"} {
        if (NR == 1) { print \$0, "Type"; next }
        if ((length(\$4) > 1 || length(\$5) > 1) && length(\$4) != length(\$5)) { \$9="3" }
        else if (\$8 == "1") { \$9="1" }
        else if (\$8 == "0/1" || \$8 == "1/0" || \$8 == "0|1" || \$8 == "1|0" ) { \$9="2" }
        else { \$9="UNKNOWN" }
        print
    }' ${vcf_file.baseName}.${method}.filtered.tmp.txt > ${vcf_file.baseName}.${method}.filtered.txt

    rm ${vcf_file.baseName}.${method}.filtered.tmp.txt

    """
}

