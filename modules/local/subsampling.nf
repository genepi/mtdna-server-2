process SUBSAMPLING {

    input:
    path bam_file
    val coverage

    output:
    path "${bam_file}", includeInputs: true, emit: subsampled_bam_ch
    
    """
    # calculate mean coverage
    samtools coverage ${bam_file} > samtools_coverage_${bam_file.baseName}.txt
    csvtk grep -t -f3 -p 16569 -C '\$' samtools_coverage_${bam_file.baseName}.txt
    mean_cov=\$(csvtk grep -t -f3 -p 16569 -C '\$' samtools_coverage_${bam_file.baseName}.txt | csvtk cut -t -f 7)
    
    # convert to integer
    mean_cov_int=\$(printf "%.0f" "\$mean_cov")
    
    # set seed to 1 (for reproducbility) and FRAC to coverage/mean_cov (e.g. 1.2)
    fraction=\$(echo "scale=4; 1+(${coverage} / \${mean_cov})" | bc)
    
    if [ \${mean_cov_int} -gt ${coverage} ]
    then
        samtools view -s \$fraction -b -o ${bam_file.baseName}.subsampled.bam ${bam_file}
        mv ${bam_file.baseName}.subsampled.bam ${bam_file}
    fi 
    """
}