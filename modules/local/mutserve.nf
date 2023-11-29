process MUTSERVE {

    publishDir "${params.output}/variants", mode: 'copy'

    input:
    path bam_files_ch
    path reference

    output:
    path("variants.vcf.gz"), emit: vcf_ch
    path("variants.txt"), emit: txt_ch

    """
    java -jar /opt/mutserve/mutserve.jar call --level 0.01 --reference ${reference} --mapQ 30 --baseQ 30 --deletions --output variants.vcf.gz --no-ansi ${bam_files_ch} --threads 1 --write-raw
    """
}