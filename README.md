ExtractorFlow pipeline
======================

Usage
-----

See *USAGE* or run `nextflow run main.nf --help`

Singularity
-----------

If you are on Linux, we recommend using the Singularity container to run ExtractorFlow

Run the command (your working directory has to be the "containers" directory):

`sudo singularity build singularity_extractorflow.sif singularity_extractorflow.def`

Then you can start the pipeline using this command line:

```
nextflow main.nf --root=/path_to_your_data/ -with-singularity singularity_extractorflow.sif -resume
```


Docker
------
If you are on MacOS or Windows, we recommend using the Docker container to run ExtractorFlow.

You can build docker image using this command (your working directory has to be the "containers" directory):

`docker build -t extractor_flow .`

```
nextflow main.nf --root=/path_to_your_data/ -with-docker extractor_flow:latest -resume
```
