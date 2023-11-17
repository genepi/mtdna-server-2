# nf-mito-calling

> Nextflow pipeline to analyze mtDNA heteroplasmy

## Requirements

- Nextflow:

```
curl -s https://get.nextflow.io | bash
```

- Docker

## Installation

Build docker image before run the pipeline:

```
docker build -t haansi/nf-mito-calling . # don't ingore the dot here
```

Test the pipeline and the created docker image with test-data:

```
nextflow run main_bam.nf
```

## Usage

```
nextflow run main_bam.nf --input input/bam/bamfile*.bam --output results
```

## Contact

- Hansi Weissensteiner (@haansi), Institute of Genetic Epidemiology, Medical University of Innsbruck


## License
MIT Licensed.
