FROM ubuntu:latest

MAINTAINER jiaxin.wang@sjtu.edu.cn

RUN apt-get update
RUN apt-get install -y --fix-missing build-essential time gfortran libblas-dev liblapack-dev gsl-bin libgsl-dev autoconf cmake
RUN apt-get install -y --fix-missing git wget vim unzip valgrind python3 python3-dev python3-pip
RUN pip3 install numpy scipy matplotlib astropy healpy jupyter pytest nbval

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
RUN wget http://downloads.sourceforge.net/project/healpix/Healpix_3.50/Healpix_3.50_2018Dec10.tar.gz && tar xzf Healpix_3.50_2018Dec10.tar.gz -C healpix --strip-components 1
WORKDIR /home/lab/healpix
# write config file on fly
RUN echo '3\n\
	gfortran\n\
	\n\
	y\n\
	-I$(F90_INCDIR) -DGFORTRAN -fno-second-underscore\n\
	-O3\n\
	gcc\n\
	-O3 -std=c99 -DgFortran\n\
	ar -rsv\n\
	libcfitsio.a\n\
	/usr/local/lib\n\
	n\n\
	0\n\
	y\n\
	y\n\
	y\n\
	0\n'\
	> hlpx_config
RUN ./configure -L < hlpx_config && make
ENV LD_LIBRARY_PATH /home/lab/healpix/lib:${LD_LIBRARY_PATH}
WORKDIR /home/lab

# GOOGLE TEST
RUN git clone https://github.com/google/googletest.git googletest
WORKDIR /home/lab/googletest
RUN mkdir /home/lab/googletest/build
WORKDIR /home/lab/googletest/build
RUN cmake -DCMAKE_INSTALL_PREFIX:PATH=/usr/local .. && make && make install && cp ../googletest/src/* /usr/local/src
WORKDIR /home/lab

# MIFEOS
RUN git clone https://github.com/gioacchinowang/mifeos.git mifeos

# CLEAN PACKAGES
RUN rm -f *.tar.gz
# SET PATH
ENV LD_LIBRARY_PATH /usr/local/lib:${LD_LIBRARY_PATH}