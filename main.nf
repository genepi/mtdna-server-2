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

include { MITOCALLING } from './workflows/mitocalling'

workflow {
    MITOCALLING ()
}

