

include { MUTSERVE } from '../modules/local/mutserve'

workflow MITOCALLING {

    println "Welcome to ${params.service.name}"

    bams_ch = Channel.fromPath(params.input)
    ref_file= file(params.reference)

    MUTSERVE(
        bams_ch,
        ref_file
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