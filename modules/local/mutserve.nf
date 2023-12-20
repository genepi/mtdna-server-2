process MUTSERVE {

    publishDir "${params.output}", mode: 'copy',pattern: '*.vcf.gz'

    input:
    path bam_files_ch
    path reference
    path excluded_samples

    output:
    path("variants.vcf.gz"), emit: vcf_ch
    path("variants.txt"), emit: txt_ch

    script:
    def avail_mem = 1024
    if (!task.memory) {
        log.info '[MUTSERVE] Available memory not known - defaulting to 1GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }    

    """
    java -Xmx${avail_mem}M -jar /opt/mutserve/mutserve.jar call \
    --level ${params.detection_limit} \
    --reference ${reference} \
    --mapQ ${params.mapQ} \
    --baseQ ${params.baseQ} \
    --deletions --output variants.vcf.gz \
    --no-ansi ${bam_files_ch} \
    --excluded-samples ${excluded_samples} \
    --threads ${task.cpus} \
    --write-raw
    """
}