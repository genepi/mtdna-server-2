# mtDNA-Server 2

[![mtdna-server-2](https://github.com/genepi/mtdna-server-2/actions/workflows/run-tests.yml/badge.svg)](https://github.com/genepi/mtdna-server-2/actions/workflows/run-tests.yml)
[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://github.com/askimed/nf-test)

mtDNA-Server 2 is based on Nextflow and accurately detects heteroplasmic and homoplasmic variants in mitochondrial genomes. It supports different variant callers and is able to call insertions, deletions and single nucleotide variants at once. Furthmore it includes (a) quality control and input validation modules, (b) haplogroup classification and contamination detected, (c) a method to estimate the individual required coverage to minimize false positives and (d) an interactive analytics dashboard. 

![image](docs/images/report.png)

## Web Service

mtDNA-Server 2 is hosted as a **free** service on our [mitoverse](https://mitoverse.i-med.ac.at/) platform.

## Documentation

Documentation can be accessed [here](https://mitoverse.readthedocs.io/mtdna-server/mtdna-server/).

![image](docs/images/workflow.png)

## Command-Line Execution with Nextflow 

### Requirements

- Docker or Singularity
- Java
- Nextflow

```
nextflow run genepi/mtdna-server-2 -r v2.1.6 -profile test,docker
```

## Publication

Weissensteiner H*, Forer L*, Fuchsberger C, Schöpf B, Kloss-Brandstätter A, Specht G, Kronenberg F, Schönherr S: mtDNA-Server: next-generation sequencing data analysis of human mitochondrial DNA in the cloud. Nucleic Acids Res. 44:W64-9, 2016. PMID: [27084948](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4987870/)

## Contact

This software was developed at the [Institute of Genetic Epidemiology](https://genepi.i-med.ac.at/), [Medical University of Innsbruck](https://i-med.ac.at/)

![](https://avatars2.githubusercontent.com/u/1942824?s=30) [Sebastian Schoenherr](mailto:sebastian.schoenherr@i-med.ac.at) ([@seppinho](https://twitter.com/seppinho))

![](https://avatars2.githubusercontent.com/u/1931865?s=30) [Hansi Weissensteiner](mailto:hansi.weissensteiner@i-med.ac.at) ([@whansi](https://twitter.com/whansi))

![](https://avatars2.githubusercontent.com/u/210220?s=30) [Lukas Forer](mailto:lukas.forer@i-med.ac.at) ([@lukfor](https://twitter.com/lukfor))

## License

MIT Licensed.
