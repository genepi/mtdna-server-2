process MUTECT2 {

    publishDir "${params.output}", mode: 'copy'

    input:
    path bam_file
    path reference
    path excluded_samples
    path fasta_index_files
    val detected_contig

    output:
    path("${bam_file.baseName}.vcf.gz"), emit: vcf_ch
    path("${bam_file.baseName}.txt"), emit: txt_ch

    """
    samtools index ${bam_file}

    gatk Mutect2 \
    --mitochondria-mode \
    -R ${reference} \
    -L ${detected_contig} \
    --min-base-quality-score ${params.mutserve.baseQ} \
    -I ${bam_file} \
    -O raw.vcf.gz
    
    gatk FilterMutectCalls \
    --mitochondria-mode \
    --max-alt-allele-count 10 \
    -R ${reference} \
    --min-allele-fraction ${params.detection_limit} \
    -V raw.vcf.gz \
    -O ${bam_file.baseName}.vcf.gz

    rm raw.vcf.gz

    echo -e "ID\tFilter\tPos\tRef\tVariant\tVariantLevel\tCoverage" > ${bam_file.baseName}.txt

    bcftools query -f '${bam_file}\t%FILTER\t%POS\t%REF\t%ALT\t[%AF\t%DP]\n' ${bam_file.baseName}.vcf.gz >> ${bam_file.baseName}.txt

    """
}