FROM continuumio/miniconda3
RUN mkdir -p /usr/local/app
WORKDIR /usr/local/app
COPY . ./
VOLUME /usr/local/app/results/aurora_test_real2
RUN apt update -y && apt install -y build-essential ffmpeg && apt clean
RUN conda env create -f envs/env_aurora1.yml
CMD ["./forecast.sh"]
