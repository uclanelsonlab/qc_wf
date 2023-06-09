# docker build -t qualimap:2.2.2a . --platform linux/amd64

FROM continuumio/miniconda3:23.3.1-0-alpine

# Install libfreetype6
#  RUN apk upgrade --update && \
#     apk add libfreetype6-dev
# Install qualimap v2.2.2a 
RUN conda install -c bioconda qualimap=2.2.2a -y