process CALCULATE_STATISTICS {

    input:
    path bam_file

    output:
    path "*summary.tab.txt", emit: stats_ch
    path "*.zip", emit: fastqc_ch


    script:
    def output_name = "${bam_file.baseName}.summary.txt"
    def output_name_tab = "${bam_file.baseName}.summary.tab.txt"

    """
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
    echo Sample\tParameter\tValue > $output_name
    echo ${bam_file.baseName}\tContig\t\${contig} >> $output_name
    echo ${bam_file.baseName}\tNumberofReads\t\${numreads} >> $output_name
    echo ${bam_file.baseName}\tCoveredBases\t\${covered_bases} >> $output_name
    echo ${bam_file.baseName}\tCoveragePercentage\t\${covered_bases_percentage} >> $output_name
    echo ${bam_file.baseName}\tMeanDepth\t\${mean_depth} >> $output_name
    echo ${bam_file.baseName}\tMeanBaseQuality\t\${mean_base_quality} >> $output_name
    echo ${bam_file.baseName}\tMeanMapQuality\t\${mean_map_quality} >> $output_name

    csvtk space2tab $output_name -o $output_name_tab

    fastqc $bam_file -o .

    """
}