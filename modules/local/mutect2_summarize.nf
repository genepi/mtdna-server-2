process MUTECT2_SUMMARIZE {

    input:
    path variants_txt

    output:
    path("variants.txt"), emit: txt_summarized_ch

    """
    csvtk concat -t ${variants_txt} -T -o variants.txt
    """
}