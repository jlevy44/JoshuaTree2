FROM continuumio/miniconda:latest

RUN pip install --upgrade pip

RUN apt-get update && apt-get install -y build-essential

RUN pip install biopython matplotlib git+git://github.com/tanghaibao/jcvi.git numpy pandas

RUN pip install --pre sourmash

RUN conda install -y -c bioconda lastz bbmap bedtools gffread gmap

ENTRYPOINT ["/bin/sh", "-c"]