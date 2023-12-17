# mtDNA-Server 2
[![mtdna-server-2](https://github.com/genepi/mtdna-server-2/actions/workflows/run-tests.yml/badge.svg)](https://github.com/genepi/mtdna-server-2/actions/workflows/run-tests.yml)
[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://github.com/askimed/nf-test)

This repository includes the mtDNA-Server-2 Workflow ported to Nextflow.

## Access mtDNA-Server 2
- mtDNA-Server 2 is hosted as a service on [mitoverse](https://mitoverse.i-med.ac.at/)
- mtDNA-Server 2 can also be installed [locally](##mtdna-server-2) as a service.  
## Requirements

- Nextflow:

```
curl -s https://get.nextflow.io | bash
```


## Run mtdna-Server 2 locally

Build docker image before run the pipeline:

```
docker build -t genepi/mtdna-server-2 . # don't ingore the dot here
nextflow run main.nf -c tests/test_single_bam.config -profile development


## Run mtdna-Server 2 as a webservice

### Requirements:

- Install Nextflow
- Docker or Singularity
- Java 14

### Installation

- Install cloudgene3: `curl -s install.cloudgene.io | bash -s 3.0.0-beta4`
- Download latest source code zip file from releases
- Install nf-impuationserver app: `./cloudgene install mtdna-server2.zip`
- Start cloudgene server: `./cloudgene server`
- Open [http://localhost:8082](http://localhost:8082)
- Login with default admin account: username `admin` and password `admin1978`
- mtDNA-Server-2 can be tested with the following [test file](https://github.com/genepi/mtdna-server-2/raw/main/tests/data/bam/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam))


## Contact
Institute of Genetic Epidemiology, Medical University of Innsbruck 

## License
MIT Licensed.
