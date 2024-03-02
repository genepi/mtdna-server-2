process SAMPLE_REPORT {

    publishDir "${params.output_reports}/samples", mode: 'copy'

    input:
    path report
    path variants
    path haplogroups
    path haplocheck
    path statistics
    path mapping
    path excluded
    tuple val(sample_filename), val(sample), file(sample_bam_file)

    output:
    file "*.html" 

    """
    echo -e "Parameter\tValue" > params.txt
    echo -e "Version\t${workflow.manifest.version}" >> params.txt
    echo -e "Job\t${params.project}" >> params.txt
    echo -e "Date\t${params.project_date}" >> params.txt  
    echo -e "Repository\t${params.service.github}" >> params.txt    
    echo -e "Variant Caller\t${params.mode}" >> params.txt
    echo -e "Detection Limit\t${params.detection_limit}" >> params.txt
    echo -e "Reference\t${params.reference}" >> params.txt
    echo -e "Base Quality\t${params.baseQ}" >> params.txt
    echo -e "Map Quality\t${params.mapQ}" >> params.txt
    echo -e "Alignment Quality\t${params.alignQ}" >> params.txt  

    #TODO: split fwd and reverse? https://bioinformatics.stackexchange.com/questions/8649/how-can-i-calculate-coverage-at-single-bases-using-a-bam-file
    echo -e "Contig\tPosition\tCoverage" > coverage.tsv
    samtools depth -a  ${sample_bam_file} >> coverage.tsv

    Rscript -e "require('rmarkdown'); render('${report}',
    params = list(
        pipeline_parameters = 'params.txt',
        variants = '${variants}',
        haplogroups = '${haplogroups}',
        haplocheck = '${haplocheck}',
        statistics = '${statistics}',
        mapping = '${mapping}',
        excluded_samples = '${excluded}',
        sample = '${sample}',
        coverage = 'coverage.tsv'
    ),
    knit_root_dir='\$PWD', output_file='\$PWD/${sample}.html')"
    """

}