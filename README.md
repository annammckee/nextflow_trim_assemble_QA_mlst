# NextFlow Trim Assemble QA and Sequence Type
This nextflow script takes paired end fastq.gz files trims and assembles them with trimmomatic and SPAdes. Then conducts QA with QUAST and characterizes bacterial isolates with MLST

#### This script trims fastq.gz files with trimmoatic, assembles them with SPAdes, runs quality assurance with QUAST, and does multi-locus sequence typing with MLST

##### Create conda environment with software
```
conda env create -f quast.yml
```

##### Activate conda environment
```
conda activate quast
```

##### Run nextflow script
```
nextflow run biol7210_nf_homework2.nf
```
