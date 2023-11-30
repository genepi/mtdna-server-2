#!/usr/bin/env nextflow

params.input = 'input/fastq/*.fastq.gz'
params.threads = 4
params.reference ="input/reference/chrM.fasta"
params.output = "output"

sequences_ch = Channel.fromPath(params.input)

ref_ch= Channel.fromPath(params.reference)
ref_ch.into { ref_mmi; ref_mutserve }

report_file = file("$baseDir/reports/report4mutserve2.rmd")


process generateMMI{
  input:
  file reference from ref_mmi

  output:
  path "${reference.baseName}.mmi" into refMMI_ch

  """
  minimap2 -k 18 -d ${reference.baseName}.mmi $reference
  """
}

process alignMap {

    cpus params.threads

    input:
    path seq from sequences_ch
    path refMMI from refMMI_ch.collect()

    output:
    file "${params.output}.sam" into samFile_ch

    """
    echo $seq
    minimap2 -a  $refMMI $seq -t ${task.cpus} > ${params.output}.sam
    """
}

process bwa_sort_index {
    publishDir "$params.output", mode: 'copy'

    cpus params.threads

    input:
    file sam_file from samFile_ch

    output:
    file "${sam_file.baseName}.bam" into sorted_bam
    file "${sam_file.baseName}.bam.bai" into sorted_bam_index

    """
    samtools sort -o ${sam_file.baseName}.bam -O BAM $sam_file -@ ${task.cpus}
    samtools index ${sam_file.baseName}.bam
    rm $sam_file
    """
}

process run_mutserve_call {
    publishDir "$params.output", mode: 'copy'

    input:
    file bam_file from sorted_bam
    file reference from ref_mutserve

    output:
    file "${params.output}.vcf" into mutserve_out_vcf
    file "${params.output}.txt" into mutserve_out_txt
    file "${params.output}_raw.txt" into mutserve_out_raw

    """
    mutserve call --output ${params.output}.vcf --reference $reference $bam_file --write-raw --threads ${task.cpus}
    """
}

process run_mutserve_annotation {
    publishDir "$params.output", mode: 'copy'

    input:
    file variants from mutserve_out_txt

    output:

    file "${params.output}_annotated.txt" into mutserve_ann_variants

    """
    mutserve annotate --input $variants --annotation /opt/mutserve/rCRS_annotation_2020-08-20.txt --output  ${params.output}_annotated.txt
    """
}


process run_haplocheck{
    publishDir "$params.output", mode: 'copy'

    input:
    file vcf_file from mutserve_out_vcf

    output:
    file "${params.output}_contaminated.txt" into haplocheck_out
    file "${params.output}_contaminated.raw.txt" into haplocheck_ext_out

    """
    haplocheck --out ${params.output}_contaminated.txt $vcf_file --raw
    """

}


process createReport {

  publishDir "$params.output", mode: 'copy'

  input:
    file variants from mutserve_ann_variants.collect()
    file raw from mutserve_out_raw.collect()
    file contamination from haplocheck_ext_out.collect()

    file report_file

  output:
    file "*.html" into report_ch

  """
  Rscript -e "require('rmarkdown'); render('${report_file}',
   params = list(
     variants_filename = '${variants}',
     raw_filename = '${raw}',
     contamination_filename = '${contamination}'
   ),
   knit_root_dir='\$PWD', output_file='\$PWD/results.html')"
  """

}


workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
