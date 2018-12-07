# This docker file will configure an environment for FLLIT

FROM ubuntu:16.04

# Basic functionalities in Ubuntu
RUN apt-get -qq update && apt-get -qq install -y \
    xorg \
    libgomp1 \
    ghostscript \
    wget \
    unzip \
    subversion && \

# Clone the compiled FLLIT repository
# svn export https://github.com/BII-wushuang/FLLIT/trunk/Compiled /FLLIT && \

# Install the MCR dependencies
# Download and install Matlab Compiler Runtime v9.01 (2016a)
    mkdir /mcr-install && \
    mkdir /opt/mcr && \
    cd /mcr-install && \
    wget http://ssd.mathworks.com/supportfiles/downloads/R2016a/deployment_files/R2016a/installers/glnxa64/MCR_R2016a_glnxa64_installer.zip && \
    cd /mcr-install && \
    unzip -q MCR_R2016a_glnxa64_installer.zip && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf mcr-install

# Configure environment variables for MCR
ENV LD_LIBRARY_PATH /opt/mcr/v901/runtime/glnxa64:/opt/mcr/v901/bin/glnxa64:/opt/mcr/v901/sys/os/glnxa64
ENV XAPPLRESDIR /opt/mcr/v901/X11/app-defaults

CMD cd /FLLIT && bash run_FLLIT.sh /opt/mcr/v901
