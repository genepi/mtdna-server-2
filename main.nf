#!/usr/bin/env nextflow
/*
========================================================================================
    genepi/nf-mito-calling
========================================================================================
    Github : https://github.com/genepi/nf-mito-calling
    Author: Sebastian Schönherr / Hansi Weissensteiner / Lukas Forer
    ---------------------------
*/

nextflow.enable.dsl = 2


/*
========================================================================================
    RUN IMPUTATIONSERVER Workflow
========================================================================================
*/

include { MTDNA_SERVER_2 } from './workflows/mtdna_server_2'

workflow {
    MTDNA_SERVER_2 ()
}

