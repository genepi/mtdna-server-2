process INPUT_VALIDATION {

    publishDir "${params.output_auxiliary}", mode: 'copy',pattern: '*s.txt'

    input:
    path bams_ch
    path statistics
    path mapping

    output:
    path("sample_statistics.txt"), emit: summarized_ch
    path("sample_mappings.txt"), emit: mapping_ch
    path("excluded_samples.txt"), emit: excluded_ch
    path("contig.txt"), emit: contig_ch
    path("*.bam"), includeInputs: true, emit: validated_files

    """
    set +e   
    csvtk concat \
        -t ${statistics} \
        -T -o sample_statistics.txt \
        --num-cpus ${task.cpus}
    
    csvtk concat \
        -t ${mapping} \
        -T -o sample_mappings.txt \
        --num-cpus ${task.cpus}
    
    java -jar /opt/mutserve/mutserve.jar stats \
        --input sample_statistics.txt \
        --mapping sample_mappings.txt \
        --detection-limit ${params.detection_limit}  \
        --reference ${params.reference}  \
        --baseQ ${params.baseQ}\
        --mapQ ${params.mapQ} \
        --alignQ ${params.alignQ} \
        --output-excluded-samples excluded_samples.txt \
        --output-contig contig.txt \
        --report validation_report.txt \
        --min-mean-depth ${params.min_mean_coverage} \
        --min-mean-base-quality ${params.min_mean_base_quality} \
        --tool ${params.mode}
    exit_code_a=\$?

    # delete excluded_samples from BAM input channel directly
    awk -v q='"' '{print "rm " q \$1 q }' excluded_samples.txt | sh

    cat validation_report.txt
    exit \$exit_code_a
    """
}