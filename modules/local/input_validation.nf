process INPUT_VALIDATION {

    publishDir "${params.output_auxiliary}", mode: 'copy',pattern: '*s.txt'

    input:
    path statistics
    path mapping

    output:
    path("sample_statistics.txt"), emit: summarized_ch
    path("sample_mappings.txt"), emit: mapping_ch
    path("excluded_samples.txt"), emit: excluded_ch
    path("contig.txt"), emit: contig_ch

    """
    csvtk concat \
        -t ${statistics} \
        -T -o sample_statistics.txt
    
    csvtk concat \
        -t ${mapping} \
        -T -o sample_mappings.txt
    
    java -jar /opt/mutserve/mutserve.jar stats \
        --input sample_statistics.txt \
        --detection-limit ${params.detection_limit}  \
        --reference ${params.reference}  \
        --baseQ ${params.baseQ}\
        --mapQ ${params.mapQ} \
        --alignQ ${params.alignQ} \
        --output-excluded-samples excluded_samples.txt \
        --output-contig contig.txt \
        --tool ${params.mode}
    """
}