# Prokaryotic annotation tools benchmarking 
This repository contains the code used in the following work:

    TODO  [Insert publication link or details here]
    
## Short Description 
This Nextflow workflow was developed to benchmark the performance of four open-source prokaryotic genome annotation tools: Prokka, Bakta, EggNOG-mapper, and PGAP. We evaluated these tools on 180,504 diverse genomes to provide guidance on tool selection based on genome characteristics and various investigated metrics. For furter details please look at the publication.

## Installation

### Dependencies

* [Nextflow](https://www.nextflow.io/index.html)
* [Docker](https://docs.docker.com/install/)
### Tools & Databases 
The tools and their respective database will be pulled automaticaly. 
  - Prokka: [Prokka 1.14.6](https://github.com/tseemann/prokka) with its corresponding database version.
  - Bakta:  [Bakta 1.7](https://github.com/oschwengers/bakta) with [database version 5.0](https://zenodo.org/records/7669534) 
  - EggNOG-mapper: [EggNOG-mapper 2.1.9](https://github.com/eggnogdb/eggnog-mapper) with its corresponding database version.
  - PGAP: [PGAP 2022-10-03.build 6384](https://github.com/ncbi/pgap) with its corresponding database version.

### Disclamer
As of last checked on July 8, 2025, the corresponding version of the PGAP database (https://s3.amazonaws.com/pgap/) is not anymore available for download.
### Execution 
To run the workflow, use the following command:
````bash
nextflow run  AggressiveHayBale/anno_bench  --csv sample_list.csv -profile local,docker
````
The sample_list.csv file should be formatted as follows:
```csv
Accession,Species,Path,Noise,Noise2,Noise3
GCA_000016605.1,Metallosphaera sedula,[path...]/GCA_000016605.1.fasta,0.005,0.01,0.02
```
Where: 
- Accession - NCBI GenBank Assembly Accession
- Species - NCBI-compliant taxonomy (genus or species name).
- Path - Absolute path to a FASTA genome file.
- Noise - Percentage of introduced frameshifts.

  ### Results
The pipeline output will be stored in the "results" folder.
