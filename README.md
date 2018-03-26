# Fog Simulation on Real Scenes for Partially Synthetic Foggy Data

Created by Christos Sakaridis at Computer Vision Lab, ETH Zurich.


### Introduction

This is the source code for the fog simulation pipeline that we present in our IJCV article titled [**Semantic Foggy Scene Understanding with Synthetic Data**][project_page], which is used to create the *Foggy Cityscapes* dataset from the original [Cityscapes][cityscapes] dataset. This pipeline involves computation of a denoised and complete depth map as well as a corresponding transmittance map, the latter being combined with the real clear-weather image to obtain the partially synthetic foggy image.


### Citation

If you use our fog simulation code in your work, please cite:
- our publication as listed on our [website][project_page],
- the [Cityscapes publication][cityscapes_citation], and
- the [SLIC superpixels publication][slic_citation].


### License

Our fog simulation code is made available for non-commercial use under the license agreement which is contained in the LICENSE.txt file.


### Contents

1. [Requirements](#requirements)
3. [Basic installation](#installation-for-running-the-demo)
4. Demo
5. Beyond the demo: fog simulation at dataset scale
6. Extra downloads: **CNN models** fine-tuned on *Foggy Cityscapes*


### Requirements

1.  **MATLAB**: required
    The code has been tested with MATLAB releases 2016b and 2017b. We therefore recommend using release 2016b or later, noting that these releases are only available for 64-bit systems. If such a configuration is not possible, using an earlier MATLAB release is the recommended (though not tested) alternative.
2.  C compiler: required for compilation of binary MEX file for SLIC superpixels
3.  Linux 64-bit: recommended
    A binary MEX file for SLIC superpixels is provided in this repository for Linux 64-bit platforms (file extension: `mexa64`). Users working on other platforms have to create the respective MEX file themselves via MATLAB (see [instructions](#installation-for-running-the-demo) below).


### Installation for running the demo


### Dataset Structure

The folder structure of *Foggy Cityscapes* follows that of Cityscapes, with minor extensions and modifications. It is summarized as follows:
```
{root}/{type}/{split}/{city}/{city}_{seq:0>6}_{frame:0>6}_{type}{ext}
```
Please refer to the README file of the Cityscapes git repository for a detailed presentation of this structure: 
https://github.com/mcordts/cityscapesScripts
In the following, we outline only those elements of the folder structure that are potentially differentiated in *Foggy Cityscapes*.

The meaning of the individual elements is:
 - `root`  the root folder of the Foggy Cityscapes dataset. We recommend that this folder coincides with the root folder of the original Cityscapes dataset.
 - `type`  the type/modality of data, e.g. `depth_stereoscopic` for completed depth maps, or `leftImg8bit` for left 8-bit foggy images.
 - `split` the split, i.e. `train`/`val`/`test`/`train_extra`.
 - `ext`   the extension of the file and optionally a suffix, e.g. `.mat` for completed depth maps.

Possible values of `type`
 - `leftImg8bit`               the left foggy images in 8-bit format. These are the foggy versions of the standard annotated Cityscapes images.
 - `depth_stereoscopic`        denoised and complete depth maps, based on the precomputed disparity maps (`disparity` type) included in Cityscapes. These depth maps are stored in MATLAB `.mat` files in `double` format and **directly encode depth values in meters**.
 - `leftImg8bit_transmittance` the transmittance maps for the left foggy images, in 8-bit format. Value `0` corresponds to transmittance `0`, while value `255` corresponds to transmittance `1`.


### Foggy Cityscapes-refined

A refined list of 550 Cityscapes images (498 `train` plus 52 `val`) that yield high-quality synthetic foggy images is provided in the file `foggy_trainval_refined_filenames.txt`. Details on the refinement criteria are given in our [publication][project_page].


### Contact

Christos Sakaridis
csakarid[at]vision.ee.ethz.ch
people.ee.ethz.ch/~csakarid/SFSU_synthetic

[project_page]: <http://people.ee.ethz.ch/~csakarid/SFSU_synthetic/>
[cityscapes]: <https://www.cityscapes-dataset.com/>
[cityscapes_citation]: <https://www.cityscapes-dataset.com/citation/>
[slic_citation]: <https://ivrl.epfl.ch/research/superpixels/>
