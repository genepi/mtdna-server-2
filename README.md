# nf-mito-calling

> Nextflow pipeline to analyze mtDNA heteroplasmy
## Requirements

- Nextflow:

```
curl -s https://get.nextflow.io | bash
```


## Installation

Build docker image before run the pipeline:

```
docker build -t genepi/nf-mitocalling . # don't ingore the dot here
```

## Usage

```
nextflow run main.nf -c tests/data/test_single_bam.config -profile development
```

## Contact

- Hansi Weissensteiner (@haansi), Institute of Genetic Epidemiology, Medical University of Innsbruck 


## License
MIT Licensed.
