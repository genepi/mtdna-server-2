process HAPLOGROUP_DETECTION {

    publishDir "${params.output_auxiliary}", mode: 'copy'

    input:
    path variants_ch

    output:
    path ("haplogroups.txt"), emit: haplogroups_ch

    """
    java -jar /opt/haplogrep/haplogrep3.jar classify --tree phylotree-fu-rcrs@1.2 --in ${variants_ch} --out haplogroups.txt --extend-report
    """
}