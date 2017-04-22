FROM continuumio/miniconda3:4.2.12

RUN apt-get update

RUN apt-get install -y \
    gcc

RUN \
    conda config --prepend channels r &&\
    conda config --prepend channels conda-forge &&\
    conda config --prepend channels bioconda

COPY conda.pkgs /tmp/

RUN conda install --yes --file /tmp/conda.pkgs

WORKDIR /tmp/mockinbird
RUN mkdir -p {mockinbird,.git} /data
ADD .git  .git/
ADD mockinbird mockinbird/
ADD setup.* versioneer.py MANIFEST.in LICENSE.txt ./
RUN python setup.py install

WORKDIR /data
ENTRYPOINT ["mockinbird"]
