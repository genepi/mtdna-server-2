

include { MUTSERVE } from '../modules/local/mutserve'
include { ANNOTATE } from '../modules/local/annotate'
include { HAPLOGROUP_DETECTION } from '../modules/local/haplogroup_detection'
include { CONTAMINATION_DETECTION } from '../modules/local/contamination_detection'
include { REPORT } from '../modules/local/report'

workflow MITOCALLING {

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
       ref_file = file("$projectDir/files/rCRS.fasta")
       annotation_file= file("$projectDir/files/rCRS_annotation.txt")
    } else {
        exit 1, "Reference " + params.reference + "not supported"
    }

    MUTSERVE(
        bams_ch.collect(),
        ref_file
    )

    ANNOTATE(
        MUTSERVE.out.txt_ch,
        ref_file,
        annotation_file
    )

    HAPLOGROUP_DETECTION(
        MUTSERVE.out.vcf_ch
    )    

    CONTAMINATION_DETECTION(
        MUTSERVE.out.vcf_ch
    )

    REPORT(
        MUTSERVE.out.txt_ch
    )        

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