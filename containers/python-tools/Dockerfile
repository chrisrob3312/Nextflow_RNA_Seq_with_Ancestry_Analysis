FROM python:3.11-slim

LABEL description="Python tools for cancer RNA-seq pipeline"

RUN pip install --no-cache-dir \
    numpy \
    pandas \
    scipy \
    scikit-learn \
    matplotlib \
    seaborn \
    pysam \
    CrossMap

WORKDIR /opt
