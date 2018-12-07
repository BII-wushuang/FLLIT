# FLLIT
## Overview of Feature Learning-based LImb segmentation and Tracking (**FLLIT**)

The **FLLIT** program is compiled on MATLAB R2016a and runs on a Linux OS (eg. Ubuntu 16.04).

The program processes input data consisting of 512 pixels x 512 pixels video image sequences of the sample animal, e.g. *Drosophila* fruit fly. The video should be in made in a single channel (grayscale), and taken at 250 frames per second or higher speeds. Currently, TIFF and PNG files are supported. The field of view should be held steady without movement throughout the video.

**FLLIT** is able to carry out the following functions:
- Identification (at a pixel level) of the legs of the sample animal via the `Segmentation' module;
- Tracking of leg tip/claw positions via the `Tracking' module;
- Produce tracking results consisting of:
	* Body centroid position
	* Body angles of rotation (relative to the y-axis)
	* Leg trajectory in arena-centered frame of reference
	* Leg trajectory in body-centred frame of reference
- Further processing of the above raw tracking data with the `Data Process' module
- Visualisation of the tracked results via the `Make Video' module.

## Getting started with FLLIT on Ubuntu
The precompiled version of **FLLIT** can be found [here](./Compiled). Sample data is provided under the [Data](./Compiled/Data) folder. This version is compiled on MATLAB R2016a in Ubuntu 16.04 and requires the corresponding [MATLAB Runtime libraries](http://ssd.mathworks.com/supportfiles/downloads/R2016a/deployment_files/R2016a/installers/glnxa64/MCR_R2016a_glnxa64_installer.zip).

To install the MATLAB runtime libraries, open a terminal in the FLLIT directory and issue the following command

```bash install MCR_R2016a.sh```.

This will take a while to download and install the MATLAB runtime libraries to the following location:

```$HOME/MCR```.

The **FLLIT** executables consist of FLLIT and run FLLIT.sh. Open a terminal in the FLLIT directory and execute FLLIT with the following command

```bash run FLLIT.sh $HOME/MCR/v901```.

## Running FLLIT on other Operating Systems
On Windows or MacOS, it will be necessary to deploy **FLLIT** in a [Docker](https://www.docker.com/) environment. Please refer to section 1.3 and 1.4 of the [readme](./Compiled/Readme.pdf) for more details.

## Potential Issues and Troubleshooting
- It might be necessary, for first time usage, to accord executable rights to FLLIT and run FLLIT.sh, which can be done with the following command

```chmod +x run FLLIT.sh FLLIT```.

Further details about **FLLIT** can be found in the [readme](./Compiled/Readme.pdf).

A video walkthrough of **FLLIT** is covered in the following link.

[![FLLIT](https://img.youtube.com/vi/g31la3oUNYk/0.jpg)](https://www.youtube.com/watch?v=g31la3oUNYk)

