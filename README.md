# FLLIT (Published in [PLOS Biology](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000346))
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
Clone/download this repository. 

The precompiled version of **FLLIT** is located in FLLIT/Compiled.  This version is compiled on MATLAB R2016a in Ubuntu 16.04 and requires the corresponding MATLAB Runtime libraries. To install the MATLAB runtime libraries, open a terminal in the FLLIT/Compiled directory and issue the following command

```bash MCR_R2016a.sh```

This will take a while to download and install the MATLAB runtime libraries to the following location:

```$HOME/MCR```

The **FLLIT** executables consist of FLLIT and run FLLIT.sh. Open a terminal in the FLLIT directory and execute FLLIT with the following command

```bash run_FLLIT.sh $HOME/MCR/v901```

#### Sample Datasets
Sample datasets may be downloaded from [GoogleDrive](https://drive.google.com/open?id=1sHYEV85wkquBi2xYVOjj7k9g8yI5CeLE). These sample data are to be placed under the FLLIT/Compiled/Data directory.

## Running FLLIT on other Operating Systems
On Windows or MacOS, it will be necessary to deploy **FLLIT** in a [Docker](https://www.docker.com/) environment. Please refer to section 1.4 and 1.5 of the [readme](./Compiled/Readme.pdf) for more details.

## Potential Issues and Troubleshooting
1) If you face disconnection issues when downloading directly as a zip, try either of the following:

Use the git command in terminal ```git clone https://github.com/BII-wushuang/FLLIT.git```

or download with [GitHub Desktop]( https://desktop.github.com/ ) (available on Windows and MacOS).

2) It might be necessary, for first time usage, to accord executable rights to **FLLIT**, which can be done with the following command

```chmod +x FLLIT```

3) For the `Make Video' module, the [ffmpeg](https://ffmpeg.org/) libraries will be required.

- For Ubuntu, install via apt-get:

```sudo apt install ffmpeg```

- For Windows, download from the following link:

[https://ffmpeg.zeranoe.com/builds/win64/static/ffmpeg-4.3-win64-static.zip](https://ffmpeg.zeranoe.com/builds/win64/static/ffmpeg-4.3-win64-static.zip)

and either extract ffmpeg.exe to this folder or add it to environment path

- For MacOS, install via brew:

```brew install ffmpeg```

4) If the images cannot be read, please rescale them to 512x512 pixels and ensure that the image format is one of the following: ".tif", ".tiff", ".bmp", ".png", ".jpg", ".jpeg".

5) If the following error message occurs during tracking ```Operands to the || and && operators must be convertible to logical scalar values.```, please follow the steps below to manually initialize tracking

- Navigate to frame 1,     check the “**Manually Initiate Tracking**” checkbox and click on the **Tracking** button. The **Tracking** button will display **Initial**.
- Click on the **Adjust Prediction** button.  Select the legs and double click on its tip position on the image to mark its coordinates. 
- After annotating all 6 legs, click the **Save and Exit** button to record the changes. **Resume** the tracking with this manually initialized data.

### Quick Start

Further details about **FLLIT** can be found in the [readme](./Compiled/Readme.pdf).

A video walkthrough of **FLLIT** is covered in the following link.

[![FLLIT](https://img.youtube.com/vi/9846ZswUURc/0.jpg)](https://www.youtube.com/watch?v=9846ZswUURc)

