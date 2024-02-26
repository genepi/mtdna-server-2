process MUTECT2 {

    input:
    path bam_file
    path reference
    path excluded_samples
    path fasta_index_files
    val detected_contig

    output:
    path("${bam_file.baseName}.txt"), emit: mutect2_txt_ch
    tuple path("${bam_file.baseName}.vcf.gz"), val('mutect2_fusion'), emit: mutect2_fusion_vcf_ch
    path("${bam_file.baseName}.vcf.gz"), emit: mutect2_vcf_ch
    path("${bam_file.baseName}.vcf.gz.tbi"), emit: mutect2_vcf_idx_ch

    """
    samtools index ${bam_file}

    gatk Mutect2 \
        -R ${reference} \
        -L ${detected_contig} \
        --min-base-quality-score ${params.baseQ} \
        -callable-depth 6 --max-reads-per-alignment-start 0 \
        -I ${bam_file} \
        -O raw.vcf.gz
    
    gatk FilterMutectCalls \
        -R ${reference} \
        --min-reads-per-strand 2 \
        -V raw.vcf.gz \
        -O ${bam_file.baseName}.vcf.gz

    bcftools norm \
        -m-any \
        -f ${reference} \
        -o ${bam_file.baseName}.norm.vcf.gz -Oz \
        ${bam_file.baseName}.vcf.gz 
    
    mv ${bam_file.baseName}.norm.vcf.gz ${bam_file.baseName}.vcf.gz

    #required for mutect2-only mode!
    echo -e "ID\tFilter\tPos\tRef\tVariant\tVariantLevel\tCoverage\tType" \
        > ${bam_file.baseName}.txt
    
    bcftools query \
        -f '${bam_file}\t%FILTER\t%POS\t%REF\t%ALT\t[%AF\t%DP]\tMUTECT2\n' \
        ${bam_file.baseName}.vcf.gz >> ${bam_file.baseName}.txt
    
    rm raw.vcf.gz
    """
}