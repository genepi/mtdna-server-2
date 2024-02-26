println "Welcome to ${params.service.name}"

requiredParams = [
    'project', 'files', 'output'
]

for (param in requiredParams) {
    if (params[param] == null) {
        exit 1, "Parameter ${param} is required."
    }
}   
   
if (params.output_reports == null || params.output_auxiliary == null ) {
    params.output_reports = params.output
    params.output_auxiliary = params.output    
}
   
include { INDEX } from '../modules/local/index'
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
include { REPORT } from '../modules/local/report'


workflow MTDNA_SERVER_2 {
 
    report_file_ch = file("$projectDir/reports/report.Rmd")
    bams_ch = Channel.fromPath(params.files, checkIfExists:true)

    if(params.reference.equals("rcrs")){
        ref_file_mutserve = file("$projectDir/files/rcrs_mutserve.fasta")
        ref_file_mutect2 = file("$projectDir/files/mt_contigs.fasta")
        annotation_file= file("$projectDir/files/rCRS_annotation.txt")
    } else {
        exit 1, "Reference " + params.reference + "not supported"
    }

    INDEX(
        ref_file_mutect2
    )

    CALCULATE_STATISTICS(
        bams_ch
    )

    INPUT_VALIDATION(
        CALCULATE_STATISTICS.out.stats_ch.collect(),
        CALCULATE_STATISTICS.out.mapping_ch.collect()
    )


    def detected_contig = INPUT_VALIDATION.out.contig_ch.text.trim()

    QUALITY_CONTROL(
        INPUT_VALIDATION.out.excluded_ch,
        CALCULATE_STATISTICS.out.fastqc_ch.collect()
    )

    haplogrep_ch = file("$projectDir/files/haplogroups.txt")
    contamination_ch = file("$projectDir/files/haplocheck.txt")

    if (params.mode == 'mutserve') {

        MUTSERVE(
            bams_ch,
            ref_file_mutserve,
            INPUT_VALIDATION.out.excluded_ch
        )

        MERGING_VARIANTS(
            MUTSERVE.out.mutserve_txt_ch.collect(),
            params.mode
        )

        variants_txt_ch = MERGING_VARIANTS.out.txt_summarized_ch
        variants_vcf_ch = MUTSERVE.out.mutserve_vcf_ch.collect()
        variants_vcf_idx_ch = MUTSERVE.out.mutserve_vcf_idx_ch.collect()
        file_count =  MUTSERVE.out.mutserve_vcf_ch.count()
       
    } 

    else if (params.mode == 'mutect2') {

        MUTECT2(
            bams_ch,
            ref_file_mutect2,
            INPUT_VALIDATION.out.excluded_ch,
            INDEX.out.fasta_index_ch,
            detected_contig
        )

        MERGING_VARIANTS(
            MUTECT2.out.mutect2_txt_ch.collect(),
            params.mode
        )

        variants_txt_ch = MERGING_VARIANTS.out.txt_summarized_ch
        variants_vcf_ch = MUTECT2.out.mutect2_vcf_ch.collect()
        variants_vcf_idx_ch = MUTECT2.out.mutect2_vcf_idx_ch.collect()
        file_count =  MUTECT2.out.mutect2_vcf_ch.count()
    }

    else if (params.mode == 'fusion') {

        MUTSERVE(
            bams_ch,
            ref_file_mutserve,
            INPUT_VALIDATION.out.excluded_ch
        )

        MUTECT2(
            bams_ch,
            ref_file_mutect2,
            INPUT_VALIDATION.out.excluded_ch,
            INDEX.out.fasta_index_ch,
            detected_contig
        )
        
        merged_ch = MUTSERVE.out.mutserve_fusion_vcf_ch.concat(MUTECT2.out.mutect2_fusion_vcf_ch)

        FILTER_VARIANTS (
            merged_ch
        )

        MERGING_VARIANTS(
            FILTER_VARIANTS.out.combined_methods_ch.collect(),
            params.mode
        )

        variants_txt_ch = MERGING_VARIANTS.out.txt_summarized_ch
        // only use mutserve calls for haplogroup and contamination detection
        variants_vcf_ch = MUTSERVE.out.mutserve_vcf_ch.collect()
        variants_vcf_idx_ch = MUTSERVE.out.mutserve_vcf_idx_ch.collect()
        file_count =  MUTSERVE.out.mutserve_vcf_ch.count()
    }

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
