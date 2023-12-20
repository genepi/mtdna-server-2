process CALCULATE_STATISTICS {

    input:
    path bam_file

    output:
    path "*summary.txt", emit: stats_ch
    path "*mapping.txt", emit: mapping_ch
    path "*.zip", emit: fastqc_ch


    script:
    def output_name = "${bam_file.baseName}.summary.txt"
    def mapping_name = "${bam_file.baseName}.mapping.txt"

    """
    ## determine mapping
    echo -e "Sample\tFilename" > $mapping_name
    samtools samples ${bam_file} >> $mapping_name

    ## calculate summary statistics
    samtools coverage ${bam_file} > samtools_coverage_${bam_file.baseName}.txt
    csvtk grep -t -f3 -p 16569 samtools_coverage_${bam_file.baseName}.txt -T -o mtdna.txt

    contig=\$(csvtk cut -t -f 1 mtdna.txt)
    numreads=\$(csvtk cut -t -f 4 mtdna.txt)
    covered_bases=\$(csvtk cut -t -f 5 mtdna.txt)
    covered_bases_percentage=\$(csvtk cut -t -f 6 mtdna.txt)
    mean_depth=\$(csvtk cut -t -f 7 mtdna.txt)
    mean_base_quality=\$(csvtk cut -t -f 8 mtdna.txt)
    mean_map_quality=\$(csvtk cut -t -f 9 mtdna.txt)
    echo -e "Sample\tParameter\tValue" > $output_name
    echo -e "${bam_file}\tContig\t\${contig}" >> $output_name
    echo -e "${bam_file}\tNumberofReads\t\${numreads}" >> $output_name
    echo -e "${bam_file}\tCoveredBases\t\${covered_bases}" >> $output_name
    echo -e "${bam_file}\tCoveragePercentage\t\${covered_bases_percentage}" >> $output_name
    echo -e "${bam_file}\tMeanDepth\t\${mean_depth}" >> $output_name
    echo -e "${bam_file}\tMeanBaseQuality\t\${mean_base_quality}" >> $output_name
    echo -e "${bam_file}\tMeanMapQuality\t\${mean_map_quality}" >> $output_name

    fastqc $bam_file -o .

    """
}