#!/bin/bash

mkdir $HOME/mcr-install
mkdir $HOME/MCR
cd $HOME/mcr-install
wget http://ssd.mathworks.com/supportfiles/downloads/R2016a/deployment_files/R2016a/installers/glnxa64/MCR_R2016a_glnxa64_installer.zip
unzip -q MCR_R2016a_glnxa64_installer.zip
./install -destinationFolder $HOME/MCR -agreeToLicense yes -mode silent
rm -rf $HOME/mcr-install
