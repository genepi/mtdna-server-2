process MUTSERVE {

    input:
    path bam_file
    path reference

    output:
    path("${bam_file.simpleName}.txt"), emit: mutserve_txt_ch
    tuple path("${bam_file.baseName}.vcf.gz"), val('mutserve_fusion'), emit: mutserve_fusion_vcf_ch
    path("${bam_file.baseName}.vcf.gz"), emit: mutserve_vcf_ch
    path("${bam_file.baseName}.vcf.gz.tbi"), emit: mutserve_vcf_idx_ch
    
    script:
    def avail_mem = 1024
    if (task.memory) {
        avail_mem = (task.memory.mega*0.8).intValue()
    }    

    """
    java -Xmx${avail_mem}M -jar /opt/mutserve/mutserve.jar \
        call \
        --level ${params.detection_limit} \
        --reference ${reference} \
        --mapQ ${params.mapQ} \
        --baseQ ${params.baseQ} \
        --output ${bam_file.baseName}.vcf.gz \
        --no-ansi \
        --write-raw \
        ${bam_file} 

    bcftools norm \
        -m-any \
        -f ${reference} \
        -o ${bam_file.baseName}.norm.vcf.gz -Oz \
        ${bam_file.baseName}.vcf.gz 
    
    mv ${bam_file.baseName}.norm.vcf.gz ${bam_file.baseName}.vcf.gz
    tabix -f ${bam_file.baseName}.vcf.gz
    """
}