process MUTSERVE {

    publishDir "${params.output}/variants", mode: 'copy'

    input:
    path bam_files_ch
    path reference

    output:
    path("variants.vcf.gz"), emit: vcf_ch
    path("variants.txt"), emit: txt_ch

    """
    java -jar /opt/mutserve/mutserve.jar call --level ${params.detection_limit} --reference ${reference} --mapQ $params.mutserve.mapQ --baseQ $params.mutserve.baseQ --deletions --output variants.vcf.gz --no-ansi ${bam_files_ch} --threads 1 --write-raw
    """
}