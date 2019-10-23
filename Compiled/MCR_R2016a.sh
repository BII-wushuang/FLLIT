#!/bin/bash
mkdir $HOME/mcr-install
mkdir $HOME/MCR
cd $HOME/mcr-install
wget -4 "http://ssd.mathworks.com/supportfiles/downloads/R2016a/deployment_files/R2016a/installers/glnxa64/MCR_R2016a_glnxa64_installer.zip" --header "Referer: ssd.mathworks.com"
unzip -q MCR_R2016a_glnxa64_installer.zip
./install -destinationFolder $HOME/MCR -agreeToLicense yes -mode silent
rm -rf $HOME/mcr-install
