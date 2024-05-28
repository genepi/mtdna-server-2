


println "Welcome to ${params.service.name} ${workflow.manifest.version}"
println "(c) ${workflow.manifest.author}"
requiredParams = [
    'project', 'files'
]

for (param in requiredParams) {
    if (params[param] == null) {
        exit 1, "Parameter ${param} is required."
    }
}   

if(params.output == null) {
    params.pubDir = "output/${params.project}"
} else {
    params.pubDir = params.output
}

if (params.output_reports == null || params.output_auxiliary == null ) {
    params.output_reports = params.pubDir
    params.output_auxiliary = params.pubDir    
}
   
include { INDEX_CREATION } from '../modules/local/index_creation'
include { CALCULATE_STATISTICS } from '../modules/local/calculate_statistics'
include { INPUT_VALIDATION } from '../modules/local/input_validation'
include { QUALITY_CONTROL } from '../modules/local/quality_control'
include { MUTSERVE } from '../modules/local/mutserve'
include { MUTECT2 } from '../modules/local/mutect2'
include { FILTER_VARIANTS } from '../modules/local/filter_variants'
include { MERGING_VARIANTS } from '../modules/local/merging_variants'
include { VCF_MERGE } from '../modules/local/vcf_merge'
include { ANNOTATE } from '../modules/local/annotate'
include { HAPLOGROUPS_CONTAMINATION } from '../modules/local/haplogroups_contamination'
include { COVERAGE_ESTIMATION } from '../modules/local/coverage_estimation'
include { REPORT } from '../modules/local/report'
include { SUBSAMPLING } from '../modules/local/subsampling'
include { SAMPLE_REPORT } from '../modules/local/sample_report'


workflow MTDNA_SERVER_2 {
 
    report_file_ch = file("$projectDir/reports/report.Rmd", checkIfExists:true)
    sample_report_file_ch = file("$projectDir/reports/sample.Rmd", checkIfExists:true)
    bams_ch = Channel.fromPath(params.files)

    if ( params.max_samples != 0 && (bams_ch.count() > params.max_samples)) {
        def report = new CloudgeneReport()
        report.error("The maximum number of allowed samples is " + params.max_samples +".")
        exit 1
    }

    if ( bams_ch.count() == 0) {
        def report = new CloudgeneReport()
        report.error("No BAM files found.")
        exit 1
    }

    if(params.reference.equals("rcrs")){
        ref_file_mutserve = file("$projectDir/files/rcrs_mutserve.fasta")
        annotation_file= file("$projectDir/files/rCRS_annotation.txt")
    } else {
        exit 1, "Reference " + params.reference + "not supported"
    }

    if(params.subsampling.equals("on") ) {
        SUBSAMPLING (
            bams_ch,
            params.subsampling_coverage
        )
        bams_ch = SUBSAMPLING.out.subsampled_bam_ch
    }

    CALCULATE_STATISTICS(
        bams_ch
    )

    INPUT_VALIDATION(
        bams_ch.collect(),
        CALCULATE_STATISTICS.out.stats_ch.collect(),
        CALCULATE_STATISTICS.out.mapping_ch.collect()
    )

    def detected_contig = INPUT_VALIDATION.out.contig_ch.text.trim()

    INDEX_CREATION(
        ref_file_mutserve,
        detected_contig
    )

    QUALITY_CONTROL(
        CALCULATE_STATISTICS.out.fastqc_ch.collect()
    )

    haplogrep_ch = file("$projectDir/files/haplogroups.txt")
    contamination_ch = file("$projectDir/files/haplocheck.txt")

    validated_files = INPUT_VALIDATION.out.validated_files.flatten()

    if (params.mode == 'mutserve') {

        MUTSERVE(
            validated_files,
            ref_file_mutserve,
            "mutserve_single"
        )

        vcf_ch = MUTSERVE.out.mutserve_ch
        file_count =  MUTSERVE.out.mutserve_ch.count()
       
    } 

    else if (params.mode == 'mutect2') {

        MUTECT2(
            validated_files,
            INDEX_CREATION.out.ref_ch,
            INDEX_CREATION.out.fasta_index_ch,
            detected_contig,
            "mutect2_single"
        )

        vcf_ch = MUTECT2.out.mutect2_ch
        file_count =  MUTECT2.out.mutect2_ch.count()
    }

    else if (params.mode == 'fusion') {

        MUTSERVE(
            validated_files,
            ref_file_mutserve,
            "mutserve_fusion"
        )

        MUTECT2(
            validated_files,
            INDEX_CREATION.out.ref_ch,
            INDEX_CREATION.out.fasta_index_ch,
            detected_contig,
            "mutect2_fusion"
        )
        
        vcf_ch = MUTSERVE.out.mutserve_ch.concat(MUTECT2.out.mutect2_ch)
        file_count =  MUTSERVE.out.mutserve_ch.count()
    }

    FILTER_VARIANTS (
        vcf_ch
    )

    MERGING_VARIANTS(
        FILTER_VARIANTS.out.combined_methods_ch.collect(),
        params.mode
    )
    
    variants_txt_ch = MERGING_VARIANTS.out.txt_summarized_ch    

    if (params.mode != 'mutect2') {

        variants_vcf_ch = vcf_ch
            .filter { it[2].contains("mutserve") }
            .collect { it[0]}

        variants_vcf_idx_ch = vcf_ch
            .filter { it[2].contains("mutserve") }
            .collect { it[1]}

        VCF_MERGE (
            variants_vcf_ch,
            variants_vcf_idx_ch,
            file_count
        )

        HAPLOGROUPS_CONTAMINATION (
            VCF_MERGE.out.vcf_merged_ch
        )
        haplogrep_ch = HAPLOGROUPS_CONTAMINATION.out.haplogroups_ch
        contamination_ch =  HAPLOGROUPS_CONTAMINATION.out.contamination_txt_ch
    }

    if (params.coverage_estimation.equals("on") && params.mode != 'mutect2') {
        COVERAGE_ESTIMATION(
            variants_txt_ch
        )
        variants_txt_ch = COVERAGE_ESTIMATION.out.variants_verified_ch
    }    
    
    ANNOTATE(
        variants_txt_ch,
        ref_file_mutserve,
        annotation_file
    )

    REPORT(
        report_file_ch,
        variants_txt_ch,
        haplogrep_ch,
        contamination_ch,
        INPUT_VALIDATION.out.summarized_ch,
        INPUT_VALIDATION.out.mapping_ch,
        INPUT_VALIDATION.out.excluded_ch
    )

    bams_filenames = bams_ch
        .map (bam_file -> tuple(bam_file.name, file(bam_file)))

    samples = INPUT_VALIDATION.out.mapping_ch
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.Filename, row.Sample) }
        .join(bams_filenames)


    SAMPLE_REPORT(
        sample_report_file_ch,
        variants_txt_ch,
        haplogrep_ch,
        contamination_ch,
        INPUT_VALIDATION.out.summarized_ch,
        INPUT_VALIDATION.out.mapping_ch,
        INPUT_VALIDATION.out.excluded_ch,
        samples
    )


}

workflow.onComplete {
  
    def report = new CloudgeneReport()
   
    //job failed
    if (!workflow.success) {
        if (params.config.send_mail){
            sendMail{
                to "${params.user.email}"
                from "${params.service.email}"
                subject "[${params.service.name}] Job ${params.project} failed."
                body "Hi ${params.user.name}, your job failed :(. Logs can be accessed at ${params.service.url}/index.html#!jobs/${params.project}"

            }
        }
        report.error("Job failed. Reason: " + workflow.errorMessage)
        return
    } else {
        report.ok("Job terminated successfully. Duration: " + workflow.duration)
    }

    //job successful
    if (params.config.send_mail){
        sendMail{
            to "${params.user.email}"
            from "${params.service.email}"
            subject "[${params.service.name}] Job ${params.project} is complete."
            body "Hi ${params.user.name}, your job completed successfully and can be accessed at ${params.service.url}/index.html#!jobs/${params.project}"
        }
        report.ok("Sent email notification to <b>${params.user.email}</b>")
    } else {

    }
}