process REPORT {

  publishDir "$params.output", mode: 'copy'

  input:
    path report
    path variants
    path haplogroups
    path haplocheck
    path statistics

  output:
    file "*.html" 

  """
  Rscript -e "require('rmarkdown'); render('${report}',
   params = list(
       variants = '${variants}',
       haplogroups = '${haplogroups}',
       haplocheck = '${haplocheck}',
       statistics = '${statistics}'
   ),
   knit_root_dir='\$PWD', output_file='\$PWD/results.html')"
  """

}