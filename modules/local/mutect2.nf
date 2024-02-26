process MUTECT2 {

    //publishDir "${params.output}", mode: 'copy', pattern: 'vcf.gz'

    input:
    path bam_file
    path reference
    path excluded_samples
    path fasta_index_files
    val detected_contig

    output:
    tuple path("${bam_file.baseName}.vcf.gz"), val('mutect2'), emit: mutect2_vcf_ch
    path("${bam_file.baseName}.mutect2.txt"), emit: txt_ch
    path("${bam_file.baseName}.mutect2.filtered.txt"), emit: combined_results

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

    rm raw.vcf.gz

    bcftools norm -m-any -f ${reference} ${bam_file.baseName}.vcf.gz -o ${bam_file.baseName}.norm.vcf.gz -Oz
    mv ${bam_file.baseName}.norm.vcf.gz ${bam_file.baseName}.vcf.gz

    echo -e "ID\tFilter\tPos\tRef\tVariant\tVariantLevel\tCoverage\tType" > ${bam_file.baseName}.mutect2.txt

    bcftools query -f '${bam_file}\t%FILTER\t%POS\t%REF\t%ALT\t[%AF\t%DP]\tINDEL\n'  ${bam_file.baseName}.vcf.gz >> ${bam_file.baseName}.mutect2.txt

    awk -F'\t' 'NR == 1 || ((length(\$4) > 1 || length(\$5) > 1) && length(\$4) != length(\$5))' ${bam_file.baseName}.mutect2.txt > ${bam_file.baseName}.mutect2.filtered.txt

    """
}