FROM continuumio/miniconda3
RUN mkdir -p /usr/local/app
WORKDIR /usr/local/app
COPY . ./
RUN apt update -y && apt install -y build-essential ffmpeg
RUN conda env create -f envs/env_aurora1.yml
CMD ["./forecast.sh"]
