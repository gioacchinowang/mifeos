FROM ubuntu:latest

MAINTAINER jiaxin.wang@sjtu.edu.cn

RUN apt-get update
RUN apt-get install -y --fix-missing build-essential gfortran libblas-dev liblapack-dev autoconf cmake
RUN apt-get install -y --fix-missing git wget vim unzip valgrind python3 python3-dev python3-pip
RUN pip3 install numpy scipy matplotlib astropy healpy jupyter

RUN mkdir /home/lab
WORKDIR /home/lab

# CFITSIO
RUN mkdir /home/lab/cfitsio
RUN wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio_latest.tar.gz && tar xzf cfitsio_latest.tar.gz -C cfitsio --strip-components 1
WORKDIR /home/lab/cfitsio
RUN ./configure --prefix=/usr/local/ && make && make install 
WORKDIR /home/lab

# HEALPIX
RUN mkdir /home/lab/healpix
RUN wget https://sourceforge.net/projects/healpix/files/Healpix_3.60/Healpix_3.60_2019Dec18.tar.gz && tar xzf Healpix_3.60_2019Dec18.tar.gz -C healpix --strip-components 1
WORKDIR /home/lab/healpix
RUN ./configure --auto=f90
RUN make
WORKDIR /home/lab

# MIFEOS
RUN git clone https://github.com/gioacchinowang/mifeos.git mifeos

# CLEAN PACKAGES
RUN rm -f *.tar.gz
# SET PATH
ENV LD_LIBRARY_PATH /usr/local/lib:/tmp/local/lib:/home/lab/healpix/lib:${LD_LIBRARY_PATH}
