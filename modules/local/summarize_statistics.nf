process SUMMARIZE_STATISTICS {

    publishDir "${params.output}/statistics", mode: 'copy'

    input:
    path statistics

    output:
    path("sample_statistics.txt"), emit: stats_summarized_ch


    """
    csvtk concat ${statistics} -T -o sample_statistics.txt
    """
}