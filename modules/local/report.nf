process REPORT {

  publishDir "${params.output_reports}", mode: 'copy'

  input:
    path report
    path variants
    path haplogroups
    path haplocheck
    path statistics

  output:
    file "*.html" 

  """
  echo -e "Parameter\tValue" > params.txt
  echo -e "Version\t${workflow.manifest.version}" >> params.txt
  echo -e "Date\t${params.project_date}" >> params.txt  
  echo -e "Repository\t${params.service.github}" >> params.txt    
  echo -e "Variant Caller\t${params.mode}" >> params.txt
  echo -e "Detection Limit\t${params.detection_limit}" >> params.txt
  echo -e "Reference\t${params.reference}" >> params.txt
  echo -e "Base Quality\t${params.variant_calling.baseQ}" >> params.txt
  echo -e "Map Quality\t${params.variant_calling.mapQ}" >> params.txt
  echo -e "Alignment Quality\t${params.variant_calling.alignQ}" >> params.txt  

  Rscript -e "require('rmarkdown'); render('${report}',
   params = list(
       pipeline_parameters = 'params.txt',
       variants = '${variants}',
       haplogroups = '${haplogroups}',
       haplocheck = '${haplocheck}',
       statistics = '${statistics}'
   ),
   knit_root_dir='\$PWD', output_file='\$PWD/report.html')"
  """

}