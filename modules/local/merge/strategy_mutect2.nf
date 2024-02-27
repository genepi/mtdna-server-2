process STRATEGY_MUTECT2 {

    //publishDir "${params.output}", mode: 'copy'

    input:
    path bam_file
    path reference
    path excluded_samples
    path fasta_index_files
    val detected_contig

    output:
    path("${bam_file.baseName}.mutect2.filtered.txt"), emit: results
    
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

     bcftools norm -m-any -f ${reference} ${bam_file.baseName}.vcf.gz -o ${bam_file.baseName}.filtered.vcf.gz -Oz

tabix ${bam_file.baseName}.filtered.vcf.gz
    gatk LeftAlignAndTrimVariants \
   -R ${reference}  \
   -V ${bam_file.baseName}.vcf.gz \
   -O ${bam_file.baseName}.filtered2.vcf.gz \
   --max-indel-length 208

    rm raw.vcf.gz

    echo -e "ID\tFilter\tPos\tRef\tVariant\tVariantLevel\tCoverage\tType" > ${bam_file.baseName}.mutect2.txt

    bcftools query -f '${bam_file}\t%FILTER\t%POS\t%REF\t%ALT\t[%AF\t%DP]\t0\n'  ${bam_file.baseName}.filtered2.vcf.gz >> ${bam_file.baseName}.mutect2.txt

    awk -F'\t' 'NR == 1 || length(\$4) > 1 || length(\$5) > 1' ${bam_file.baseName}.mutect2.txt > ${bam_file.baseName}.mutect2.filtered.txt


    """
}