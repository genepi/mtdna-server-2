process MUTSERVE {

    input:
    path bam_file
    path reference


    """
    mutserve call --level 0.01 --reference ${reference} --mapQ 30 --baseQ 30 --deletions --output ${bam_file.baseName}.vcf.gz --no-ansi ${bam_file} --threads 1 --write-raw
    """
}