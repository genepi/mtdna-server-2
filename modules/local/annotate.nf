process ANNOTATE {

    publishDir "${params.output}", mode: 'copy'

    input:
    path variants_ch
    path reference
    path annotation

    output:
    path ("${variants_ch.baseName}_ann.txt"), emit: variants_ann_ch

    """
    java -jar /opt/mutserve/mutserve.jar annotate --input ${variants_ch} --output ${variants_ch.baseName}_ann.txt --annotation ${annotation}
    """
}
