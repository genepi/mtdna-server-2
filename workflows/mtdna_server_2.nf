include { INDEX } from '../modules/local/index'
include { CALCULATE_STATISTICS } from '../modules/local/calculate_statistics'
include { INPUT_VALIDATION } from '../modules/local/input_validation'
include { QUALITY_CONTROL } from '../modules/local/quality_control'
include { MUTSERVE } from '../modules/local/mutserve'
include { MUTECT2 } from '../modules/local/mutect2'
include { MUTECT2_SUMMARIZE } from '../modules/local/mutect2_summarize'
include { ANNOTATE } from '../modules/local/annotate'
include { HAPLOGROUP_DETECTION } from '../modules/local/haplogroup_detection'
include { CONTAMINATION_DETECTION } from '../modules/local/contamination_detection'
include { REPORT } from '../modules/local/report'

workflow MTDNA_SERVER_2 {
 
    println "Welcome to ${params.service.name}"

    requiredParams = [
        'project', 'files', 'output'
    ]

    for (param in requiredParams) {
        if (params[param] == null) {
            exit 1, "Parameter ${param} is required."
        }
    }

    bams_ch = Channel.fromPath(params.files)

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
        CALCULATE_STATISTICS.out.stats_ch.collect()
    )


    def detected_contig = INPUT_VALIDATION.out.contig_ch.text.trim()

    QUALITY_CONTROL(
        INPUT_VALIDATION.out.excluded_ch,
        CALCULATE_STATISTICS.out.fastqc_ch.collect()
    )

    if (params.mode == 'mutect2') {

        MUTECT2(
        bams_ch,
        ref_file_mutect2,
        INPUT_VALIDATION.out.excluded_ch,
        INDEX.out.fasta_index_ch,
        detected_contig
        )

        MUTECT2_SUMMARIZE(
        MUTECT2.out.txt_ch.collect()
        )

        variants_txt_ch = MUTECT2_SUMMARIZE.out.txt_summarized_ch
    }

    else {

        MUTSERVE(
            bams_ch.collect(),
            ref_file_mutserve,
            INPUT_VALIDATION.out.excluded_ch
        )

        variants_txt_ch = MUTSERVE.out.txt_ch

        HAPLOGROUP_DETECTION(
            MUTSERVE.out.vcf_ch
        )    

        CONTAMINATION_DETECTION(
            MUTSERVE.out.vcf_ch
        )
    } 

    ANNOTATE(
        variants_txt_ch,
        ref_file_mutserve,
        annotation_file
    )

    //TODO REPORT!

}

workflow.onComplete {
    //TODO: use templates
    //TODO: move in EmailHelper class
    //see https://www.nextflow.io/docs/latest/mail.html for configuration etc...
   
    def report = new CloudgeneReport()
   
    //job failed
    if (!workflow.success) {
        if (params.config.send_mail){
            sendMail{
                to "${params.user.email}"
                subject "[${params.service.name}] Job ${params.project} failed."
                body "Hi ${params.user.name}, the job failed :("
            }
        }
        report.error("Imputation failed.")
        return
    }

    //job successful
    if (params.config.send_mail){
        sendMail{
            to "${params.user.email}"
            subject "[${params.service.name}] Job ${params.project} is complete."
        }
        report.ok("Sent email with password to <b>${params.user.email}</b>")
    } else {

    }
}