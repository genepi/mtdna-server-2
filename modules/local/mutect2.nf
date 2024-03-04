process MUTECT2 {

    input:
    path bam_file
    path reference
    path fasta_index_files
    val detected_contig
    val method

    output:
    tuple path("${bam_file.baseName}.vcf.gz"), path("${bam_file.baseName}.vcf.gz.tbi"), val(method), emit: mutect2_ch

    script:
    def avail_mem = 1024
    if (task.memory) {
        avail_mem = (task.memory.mega*0.8).intValue()
    }    

    """
    samtools index ${bam_file}

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \
        Mutect2 \
        -R ${reference} \
        -L ${detected_contig} \
        --min-base-quality-score ${params.baseQ} \
        -callable-depth 6 \
        --native-pair-hmm-threads ${task.cpus} \
        --max-reads-per-alignment-start 0 \
        --tmp-dir . \
        -I ${bam_file} \
        -O raw.vcf.gz
    
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \
        FilterMutectCalls \
        -R ${reference} \
        --min-reads-per-strand 2 \
        -V raw.vcf.gz \
        --tmp-dir \$PWD \
        -O ${bam_file.baseName}.vcf.gz

    bcftools norm \
        -m-any \
        -f ${reference} \
        -o ${bam_file.baseName}.norm.vcf.gz -Oz \
        ${bam_file.baseName}.vcf.gz 
    
    mv ${bam_file.baseName}.norm.vcf.gz ${bam_file.baseName}.vcf.gz
    tabix -f ${bam_file.baseName}.vcf.gz

    rm raw.vcf.gz
    """
}