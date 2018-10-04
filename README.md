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

## Getting started with FLLIT
The precompiled version of **FLLIT** can be found [here](./Compiled). Sample data is provided under the [Data](./Compiled/Data) folder. This version is compiled on MATLAB R2016a in Ubuntu 16.04 and requires the corresponding [MATLAB Runtime libraries](http://ssd.mathworks.com/supportfiles/downloads/R2016a/deployment_files/R2016a/installers/glnxa64/MCR_R2016a_glnxa64_installer.zip).

Extract the MATLAB runtime libraries and open a terminal in the extracted folder. Initiate installation of the MATLAB runtime libraries by issuing the following command

```sudo bash install```.

By default, the library is installed to the following location:

```/usr/local/MATLAB/MATLAB Runtime/v901```.

Please take note of the **MATLAB Runtime directory** as this will be needed to execute FLLIT.

The **FLLIT** executables consist of FLLIT and run FLLIT.sh. Open a terminal in the FLLIT directory and execute FLLIT with the following command

```bash run FLLIT.sh <MATLAB Runtime directory>```.

Here the **MATLAB Runtime directory** is the installation directory of the MATLAB Runtime libraries. For example if the default installation directory is chosen, the FLLIT execution command will be

```bash run FLLIT.sh /usr/local/MATLAB/MATLAB Runtime/v901```.

## Potential Issues and Troubleshooting
- It might be necessary, for first time usage, to accord executable rights to FLLIT and run FLLIT.sh, which can be done with the following command

```chmod +x run FLLIT.sh FLLIT```.
- Be sure to have the correct location of the **MATLAB Runtime directory** when executing FLLIT.
- ```sudo``` requires administrator access. If there is no administrator access, please install the MATLAB runtime libraries without ```sudo```.

```bash ./MCR R2016a glnxa64 installer/install```.

Note that the default installation directory cannot be chosen in this case.
- In some cases, it might be required to execute FLLIT with ```sudo```

```sudo bash run FLLIT.sh <MATLAB Runtime directory>```.

Further details about **FLLIT** can be found in the [readme](./Compiled/Readme.pdf).

A video walkthrough of **FLLIT** is covered in the following link.

[![FLLIT](https://img.youtube.com/vi/g31la3oUNYk/0.jpg)](https://www.youtube.com/watch?v=g31la3oUNYk)

