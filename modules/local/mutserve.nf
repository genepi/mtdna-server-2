process MUTSERVE {

    //publishDir "${params.output}", mode: 'copy'

    input:
    path bam_file
    path reference
    path excluded_samples

    output:
    path("${bam_file.baseName}.mutserve.filtered.txt"), emit: combined_results
    tuple path("${bam_file.baseName}.vcf.gz"), val('mutserve'), emit: mutserve_vcf_ch
    tuple path("${bam_file.baseName}.vcf.gz"), path("${bam_file.baseName}.vcf.gz.tbi"), emit: mutserve_vcf_only_ch
    script:
    def avail_mem = 1024
    if (!task.memory) {
        //log.info '[MUTSERVE] Available memory not known - defaulting to 1GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }    

    """
    java -Xmx${avail_mem}M -jar /opt/mutserve/mutserve.jar call \
    --level ${params.detection_limit} \
    --reference ${reference} \
    --mapQ ${params.mapQ} \
    --baseQ ${params.baseQ} \
    --output ${bam_file.baseName}.vcf.gz \
    --no-ansi \
    --excluded-samples ${excluded_samples} \
    --write-raw \
    ${bam_file} 

    tabix ${bam_file.baseName}.vcf.gz

    echo -e "ID\tFilter\tPos\tRef\tVariant\tVariantLevel\tCoverage\tType" > ${bam_file.baseName}.mutserve.txt
    bcftools query -f '${bam_file}\t%FILTER\t%POS\t%REF\t%ALT\t[%AF\t%DP\t]SNV\n' ${bam_file.baseName}.vcf.gz >> ${bam_file.baseName}.mutserve.txt
    awk -F'\t' 'NR == 1 || (length(\$4) == 1 && length(\$5) == 1)' ${bam_file.baseName}.mutserve.txt > ${bam_file.baseName}.mutserve.filtered.txt
   
    """
}