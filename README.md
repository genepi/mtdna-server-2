# mtDNA-Server 2
[![mtdna-server-2](https://github.com/genepi/mtdna-server-2/actions/workflows/run-tests.yml/badge.svg)](https://github.com/genepi/mtdna-server-2/actions/workflows/run-tests.yml)
[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://github.com/askimed/nf-test)

Welcome to mtdna-Server 2. This repository includes the mtDNA-Server-2 Workflow ported to Nextflow. 

## Access
mtDNA-Server 2 is hosted as a service on our [mitoverse](https://mitoverse.i-med.ac.at/) platform.


## Documentation
Documentation can be accessed [here](https://mitoverse.readthedocs.io). 

## Local Installation

### Run mtDNA-Server 2 workflow

#### Requirements

- Nextflow:
```
curl -s https://get.nextflow.io | bash
```

#### Execution

```
docker build -t genepi/mtdna-server-2 . # don't ingore the dot here
nextflow run main.nf -c tests/test_single_bam.config -profile development
```

### Run mtDNA-Server 2 as a webservice

#### Requirements:

- Install Nextflow
- Docker or Singularity
- Java 14

#### Installation

- Install cloudgene3: `curl -s install.cloudgene.io | bash -s 3.0.0-beta4`
- Download latest source code zip file from releases
- Install mtDNA-Server 2 workflow: `./cloudgene install https://github.com/genepi/mtdna-server-2/archive/refs/tags/v2.0.0-beta2.zip`
- Start cloudgene server: `./cloudgene server`
- Open [http://localhost:8082](http://localhost:8082)
- Login with default admin account: username `admin` and password `admin1978`
- mtDNA-Server-2 can be tested with the following [test file](https://github.com/genepi/mtdna-server-2/raw/main/tests/data/bam/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam))

## Publication
Weissensteiner H*, Forer L*, Fuchsberger C, Schöpf B, Kloss-Brandstätter A, Specht G, Kronenberg F, Schönherr S: mtDNA-Server: next-generation sequencing data analysis of human mitochondrial DNA in the cloud. Nucleic Acids Res. 44:W64-9, 2016. PMID: [27084948](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4987870/)

## Contact
Institute of Genetic Epidemiology, Medical University of Innsbruck 

## License
MIT Licensed.
