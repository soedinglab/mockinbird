FROM continuumio/miniconda3:4.2.12

RUN \
    conda config --prepend channels r &&\
    conda config --prepend channels conda-forge &&\
    conda config --prepend channels bioconda

COPY conda.pkgs /tmp/

RUN conda install --yes --file /tmp/conda.pkgs

RUN mkdir -p /tmp/mockinbird /data
COPY mockinbird setup.py /tmp/mockinbird/
RUN pip install /tmp/mockinbird

WORKDIR /data
ENTRYPOINT ["mockinbird"]
