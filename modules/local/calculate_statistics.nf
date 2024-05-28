process CALCULATE_STATISTICS {

    input:
    path bam_file

    output:
    path "*summary.txt", emit: stats_ch
    path "*mapping.txt", emit: mapping_ch
    path "*.zip", emit: fastqc_ch
    path("*.bam"), includeInputs: true, emit: fixed_file



    script:
    def output_name = "${bam_file.baseName}.summary.txt"
    def mapping_name = "${bam_file.baseName}.mapping.txt"
 
    def avail_mem = 1024
    if (task.memory) {
        avail_mem = (task.memory.mega*0.8).intValue()
    }    
 
    """
    ## Create Mapping File
    echo -e "Sample\tFilename" > $mapping_name

    SAMPLE_NAME=\$(echo "\$(samtools samples ${bam_file})" | awk '{print \$1}')

    ## ADD RG-Tag to work with mutect2
    if [ "\$SAMPLE_NAME" == "." ]
    then
        samtools addreplacerg -r '@RG\tID:${bam_file.baseName}\tSM:${bam_file.baseName}' ${bam_file}  -o ${bam_file.baseName}_tmp.bam
        mv ${bam_file.baseName}_tmp.bam ${bam_file}
    fi

    echo "\$(samtools samples ${bam_file})" >> $mapping_name

    ## Calculate summary statistics
    samtools coverage ${bam_file} > samtools_coverage_${bam_file.baseName}.txt
    csvtk grep -t -f3 -p 16569 -C '\$' samtools_coverage_${bam_file.baseName}.txt -T -o mtdna.txt --num-cpus ${task.cpus} 
        
    contig=\$(csvtk cut -t -f 1 mtdna.txt --num-cpus ${task.cpus})
    numreads=\$(csvtk cut -t -f 4 mtdna.txt --num-cpus ${task.cpus})
    covered_bases=\$(csvtk cut -t -f 5 mtdna.txt --num-cpus ${task.cpus})
    covered_bases_percentage=\$(csvtk cut -t -f 6 mtdna.txt --num-cpus ${task.cpus})
    mean_depth=\$(csvtk cut -t -f 7 mtdna.txt --num-cpus ${task.cpus})
    mean_base_quality=\$(csvtk cut -t -f 8 mtdna.txt --num-cpus ${task.cpus})
    mean_map_quality=\$(csvtk cut -t -f 9 mtdna.txt --num-cpus ${task.cpus})
    readgroup=\$(samtools view -H ${bam_file} | csvtk grep -I -H -r -p "^@RG" --num-cpus ${task.cpus} | sed 's/\t/,/g' | head -n 1)
    
    echo -e "Sample\tParameter\tValue" > $output_name
    echo -e "${bam_file}\tContig\t\${contig}" >> $output_name
    echo -e "${bam_file}\tNumberofReads\t\${numreads}" >> $output_name
    echo -e "${bam_file}\tCoveredBases\t\${covered_bases}" >> $output_name
    echo -e "${bam_file}\tCoveragePercentage\t\${covered_bases_percentage}" >> $output_name
    echo -e "${bam_file}\tMeanDepth\t\${mean_depth}" >> $output_name
    echo -e "${bam_file}\tMeanBaseQuality\t\${mean_base_quality}" >> $output_name
    echo -e "${bam_file}\tMeanMapQuality\t\${mean_map_quality}" >> $output_name
    echo -e "${bam_file}\tRG\t\${readgroup}" >> $output_name

    fastqc --threads ${task.cpus} --memory ${avail_mem} $bam_file -o .

    """
}