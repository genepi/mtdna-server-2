# mtDNA-Server 2

[![mtdna-server-2](https://github.com/genepi/mtdna-server-2/actions/workflows/run-tests.yml/badge.svg)](https://github.com/genepi/mtdna-server-2/actions/workflows/run-tests.yml)
[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://github.com/askimed/nf-test)

mtDNA-Server 2 is a Nextflow DSL2 pipeline to accurately detect heteroplasmic and homoplasmic variants in mitochondrial (mtDNA) genomes, details can be found in our [NAR publication](https://doi.org/10.1093/nar/gkae296). 

![image](docs/images/report.png)

### Documentation
The documentation can be accessed [here](https://mitoverse.readthedocs.io/mtdna-server/mtdna-server/). 

### Running Pipeline Locally
```
git clone github.com/genepi/mtdna-server-2
nextflow run main.nf -c tests/test_single_bam.config -profile docker
````

You can find all additional DSL2 parameters and configuration options [here](https://mitoverse.readthedocs.io/mtdna-server/mtdna-server/#nextflow-dsl2-pipeline).

## Citation
Weissensteiner H*, Forer L*, Kronenberg F, Schönherr S. [mtDNA-Server 2: advancing mitochondrial DNA analysis through highly parallelized data processing and interactive analytics](https://doi.org/10.1093/nar/gkae296). Nucleic Acids Res. 2024 May 6:gkae296. doi: 10.1093/nar/gkae296. Epub ahead of print. PMID: 38709886.

## Version History

Release [v2.1.15](../../releases/tag/v2.1.15) - Add option for min mean coverage

Release [v2.1.14](../../releases/tag/v2.1.14) - Load resource conf

Release [v2.1.13](../../releases/tag/v2.1.13) - Update to latest Haplogrep3 

Release [v2.1.12](../../releases/tag/v2.1.12) - Install Haplogrep3 tree directly.

Release [v2.1.11](../../releases/tag/v2.1.11) - Improve QC command, update to latest mutserve v2.0.1.

Release [v2.1.10](../../releases/tag/v2.1.10) - Improved mutect2 support: create missing RG tags, write inidividual reference sequence on the fly, support complex ref tags.


