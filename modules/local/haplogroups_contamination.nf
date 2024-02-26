process HAPLOGROUPS_CONTAMINATION {

    publishDir "${params.output_auxiliary}", mode: 'copy'

    input:
    path(merged_vcf_ch)

    output:
    path ("haplogroups.txt"), emit: haplogroups_ch
    path ("haplocheck.html"), emit: contamination_report_ch
    path ("haplocheck.txt"), emit: contamination_txt_ch

    """
    java -jar /opt/haplogrep/haplogrep3.jar \
        classify \
        --tree phylotree-fu-rcrs@1.2 \
        --in ${merged_vcf_ch} \
        --out haplogroups.txt \
        --extend-report

    java -jar /opt/haplocheck/haplocheck.jar \
        --out haplocheck.txt \
        --raw \
        ${merged_vcf_ch}
    """
}