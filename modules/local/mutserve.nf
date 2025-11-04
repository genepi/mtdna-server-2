process MUTSERVE {

    input:
    path bam_file
    path reference
    val method

    output:
    tuple path("${bam_file.baseName}.vcf.gz"), path("${bam_file.baseName}.vcf.gz.tbi"), val(method), emit: mutserve_ch
    path("*_raw.txt", optional: true)
    
    publishDir "${params.pubDir}/mutserve_raw_files", mode: 'copy', pattern: '*_raw.txt'

    script:
    def avail_mem = 1024
    if (task.memory) {
        avail_mem = (task.memory.mega*0.8).intValue()
    }    
    def raw_option = params.mutserve_write_raw ? "--write-raw" : ""


    """
    #todo: check used mutserve strand-bias with default parameter
    samtools index ${bam_file} 
    java -Xmx${avail_mem}M -jar /opt/mutserve/mutserve.jar \
        call \
        --level ${params.detection_limit} \
        --reference ${reference} \
        --mapQ ${params.mapQ} \
        --baseQ ${params.baseQ} \
        --output ${bam_file.baseName}.vcf.gz \
        --no-ansi \
        --strand-bias 1.6 \
        $raw_option \
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