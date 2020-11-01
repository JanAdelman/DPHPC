#RUN FROM WITHIN REPO
FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

# Install base
RUN apt-get update
RUN apt-get install -y build-essential mpich

#Install google-benchmark 
RUN apt-get install -y libbenchmark-dev

RUN mkdir HPC
ADD src/ /HPC/

# docker build -t hpc .
# sudo docker run -it --entrypoint /bin/bash <ID>

#mpicxx naive_suffixarray.cpp -o main -lbenchmark
#mpiexec -np 1 ./main