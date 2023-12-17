process MUTSERVE {

    publishDir "${params.output}", mode: 'copy'

    input:
    path bam_files_ch
    path reference
    path excluded_samples

    output:
    path("variants.vcf.gz"), emit: vcf_ch
    path("variants.txt"), emit: txt_ch

    """
    java -jar /opt/mutserve/mutserve.jar call \
    --level ${params.detection_limit} \
    --reference ${reference} \
    --mapQ $params.mutserve.mapQ \
    --baseQ $params.mutserve.baseQ \
    --deletions --output variants.vcf.gz \
    --no-ansi ${bam_files_ch} \
    --excluded-samples ${excluded_samples} \
    --threads 1 --write-raw
    """
}