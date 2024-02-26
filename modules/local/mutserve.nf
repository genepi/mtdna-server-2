process MUTSERVE {

    //publishDir "${params.output}", mode: 'copy'

    input:
    path bam_file
    path reference
    path excluded_samples

    output:
    path("${bam_file.baseName}.txt"), emit: mutserve_txt_ch
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
        --excluded-samples ${excluded_samples} \
        --write-raw \
        ${bam_file} 

    bcftools norm -m-any -f ${reference} ${bam_file.baseName}.vcf.gz -o ${bam_file.baseName}.norm.vcf.gz -Oz
    mv ${bam_file.baseName}.norm.vcf.gz ${bam_file.baseName}.vcf.gz
    tabix ${bam_file.baseName}.vcf.gz
    """
}