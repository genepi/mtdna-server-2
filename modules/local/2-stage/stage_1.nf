process STAGE_1 {

    publishDir "${params.output}", mode: 'copy',pattern: '*.vcf.gz'

    input:
    path bam_file
    path reference
    path reference_mutserve
    path fasta_index_files
    val detected_contig

    output:
    path("your_output.fasta"), emit: consensus_ch

    """

    samtools index ${bam_file}

    gatk Mutect2 \
    --mitochondria-mode \
    -R ${reference} \
    -L ${detected_contig} \
    --min-base-quality-score ${params.baseQ} \
    -I ${bam_file} \
    -O raw.vcf.gz
    
    gatk FilterMutectCalls \
    --mitochondria-mode \
    --max-alt-allele-count 10 \
    -R ${reference} \
    --min-allele-fraction ${params.detection_limit} \
    -V raw.vcf.gz \
    -O ${bam_file.baseName}.vcf.gz

    bcftools consensus -H A -f ${reference} ${bam_file.baseName}.vcf.gz > consensus.fa

    samtools coverage ${bam_file} > samtools_coverage_${bam_file.baseName}.txt
    old_tag_line=\$(csvtk -H grep -t -f3 -p 16569 -C '\$' \
    samtools_coverage_${bam_file.baseName}.txt | cut -f 1)

    samtools faidx consensus.fa
    samtools faidx consensus.fa \$old_tag_line > seb.txt
    sed "s|\$old_tag_line|chrM|" seb.txt > your_output.fasta

    """
}