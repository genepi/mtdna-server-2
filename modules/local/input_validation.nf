process INPUT_VALIDATION {

    publishDir "${params.output}/statistics", mode: 'copy'

    input:
    path statistics

    output:
    path("sample_statistics.txt"), emit: summarized_ch
    path("excluded_samples.txt"), emit: excluded_ch
    path("contig.txt"), emit: contig_ch

    """
    csvtk concat -t ${statistics} -T -o sample_statistics.txt
    java -jar /opt/mutserve/mutserve.jar stats \
    --input sample_statistics.txt \
    --detection-limit ${params.detection_limit}  \
    --reference ${params.reference}  \
    --baseQ ${params.variant_calling.baseQ}\
    --mapQ ${params.variant_calling.mapQ} \
    --alignQ ${params.variant_calling.alignQ} \
    --output-excluded-samples excluded_samples.txt \
    --output-contig contig.txt \
    --tool ${params.mode}
    """
}