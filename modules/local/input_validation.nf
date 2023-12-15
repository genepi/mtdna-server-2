process INPUT_VALIDATION {

    publishDir "${params.output}/statistics", mode: 'copy'

    input:
    path statistics

    output:
    path("sample_statistics.txt"), emit: summarized_ch
    path("excluded_samples.txt"), emit: excluded_ch


    """
    csvtk concat -t ${statistics} -T -o sample_statistics.txt
    java -jar /opt/mutserve/mutserve.jar stats \
    --input sample_statistics.txt \
    --detection-limit $params.detection_limit  \
    --reference $params.reference  \
    --baseQ $params.mutserve.baseQ \
    --mapQ $params.mutserve.mapQ \
    --alignQ $params.mutserve.alignQ \
    --output excluded_samples.txt
    """
}