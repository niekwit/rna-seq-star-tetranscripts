# Dockerfiles for each version (v0.4.0 onwards)

This directory contains Dockerfiles associated with each release. 

The Docker image derived from this file contains all Conda environments for each rule, i.e. the whole workflow is run in one image.

These images are shared via [Docker Hub](https://hub.docker.com/repository/docker/niekwit/rna-seq-star-tetranscripts/general) and are generated as follows (from directory with workflow code):

```shell
$ snakemake --containerize > Dockerfile
$ sudo docker build -t niekwit/rna-seq-star-tetranscripts:{VERSION} .
$ sudo docker login
$ sudo docker push niekwit/rna-seq-star-tetranscripts:{VERSION}
```