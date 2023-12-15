process SUMMARIZE_STATISTICS {

    publishDir "${params.output}/statistics", mode: 'copy'

    input:
    path statistics

    output:
    path("sample_statistics.txt"), emit: summarized_ch
    path("excluded_samples.txt"), emit: excluded_ch


    """
    csvtk concat -t ${statistics} -T -o sample_statistics.txt
    java -jar /opt/mutserve/mutserve.jar stats --input sample_statistics.txt --baseQ $params.mutserve.baseQ --mapQ $params.mutserve.mapQ --output excluded_samples.txt
    """
}