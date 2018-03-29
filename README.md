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
4. [Demo](#demo)
5. [Beyond the demo: generating *Foggy Cityscapes*](#beyond-the-demo-running-fog-simulation-to-generate-foggy-cityscapes)
6. Extra downloads: **CNN models** fine-tuned on *Foggy Cityscapes*


### Requirements

1.  **MATLAB**  
    The code has been developed and tested in MATLAB releases 2014a, 2016b and 2017b. We therefore recommend using release 2014a or later. If such a configuration is not possible, using an earlier MATLAB release is the recommended (though not tested) alternative.
2.  **C compiler**  
    Users will have to build a binary MEX file for SLIC superpixels themselves in MATLAB (see [instructions](#installation-for-running-the-demo) below), which requires a MATLAB-supported C compiler.


### Installation for running the demo

Our fog simulation pipeline makes use of [SLIC superpixels][slic_citation], for which the algorithm is implemented in the form of a [source MEX file](source/external/SLIC_mex/slicmex.c) for usage in MATLAB. In short, a C compiler has to be configured with MATLAB and then used to build the binary MEX file for SLIC which is called in the pipeline. This [MATLAB doc page](https://www.mathworks.com/help/matlab/matlab_external/what-you-need-to-build-mex-files.html) can serve as a reference.

Steps:
1. Clone this repository with
   ```
   git clone https://github.com/sakaridis/fog_simulation-SFSU_synthetic.git
   ```
2. Make sure a MATLAB-supported C compiler is installed in your system, e.g. gcc for Linux.
3. Open MATLAB and type in the Command Window
   ```
   mex -setup
   ```
   which will guide you through the process of configuring MATLAB's `mex` command with the compiler from step 2. Once this setup has been completed successfully, issuing the above command for a second time should generate a message similar to
   ```
   MEX configured to use 'gcc' for C language compilation.
   ```
   For further details, please consult the [mex command documentation](https://www.mathworks.com/help/matlab/ref/mex.html).
4. Assign the MATLAB variable `FOG_SIMULATION_ROOT` with the path to the directory into which you have cloned this repository, by issuing in the Command Window something like
   ```
   FOG_SIMULATION_ROOT = '/home/some_user/fog_simulation-SFSU_synthetic';
   ```
5. Change MATLAB's current folder to `FOG_SIMULATION_ROOT`, e.g. with
   ```
   cd(FOG_SIMULATION_ROOT);
   ```
6. Build the binary MEX file for SLIC from the respective C source file with
   ```
   cd(fullfile('source', 'external', 'SLIC_mex'));
   mex slicmex.c;
   ```
   If the build is successful, it will generate a message similar to
   ```
   Building with 'gcc'.
   MEX completed successfully.
   ```
   and a binary MEX file `slicmex.<ext>` will be created in the [SLIC source code directory](source/external/SLIC_mex), where extension `<ext>` depends on your system (see [MATLAB docs](https://www.mathworks.com/help/matlab/matlab_external/build-an-executable-mex-file.html) for details).


### Demo

After completing the [basic installation](#installation-for-running-the-demo), you will be all set to run the demo.

To run the demo in MATLAB, type in the Command Window
```
cd(fullfile(FOG_SIMULATION_ROOT, 'source'));
Demo_fog_simulation_Cityscapes;
```

The demo runs our fog simulation on an example clear-weather image from Cityscapes and writes the results (synthesized foggy image, estimated transmittance map and depth map) under the directory `output/demos/`. The running time for a single image is around 2-3 minutes on an Intel Core i7 machine with 16 GB RAM.

We have tested the demo:
- on Linux 64-bit with gcc and MATLAB releases 2016b and 2017b
- on Windows 64-bit with Microsoft Visual C++ 2012 and MATLAB release 2014a

The results of the demo may differ for MATLAB release 2014a or earlier compared to 2014b or newer, although this difference is generally negligible. The reason is that for the required conversion of the input image to the CIE L\*a\*b\* color space, the recommended MATLAB function `rgb2lab` was only introduced in release 2014b, and the alternative implementation for earlier releases produces slightly different values for the output CIE L\*a\*b\* image. We note that for generating our *Foggy Cityscapes* dataset, we have used `rgb2lab`.


### Beyond the demo: running fog simulation to generate Foggy Cityscapes

The *Foggy Cityscapes* dataset is directly available for download at our dedicated [website][project_page] and at the [Cityscapes website][cityscapes]. However, for completeness we include in this repository some example code which can serve as a basis for users to reproduce the full-scale fog simulation experiments on Cityscapes for generating *Foggy Cityscapes*. Prior to running these experiments, the [basic installation](#installation-for-running-the-demo) must be performed and a list of Cityscapes packages must be downloaded.

1. Download from the [Cityscapes website](https://www.cityscapes-dataset.com/downloads/) the following packages:
   - `leftImg8bit_trainvaltest.zip`
   - `rightImg8bit_trainvaltest.zip`
   - `disparity_trainvaltest.zip`
   - `camera_trainvaltest.zip`
2. After unzipping these packages, you have to ensure that the extracted files obey a certain directory structure. We will denote the root directory for Cityscapes into which you have put these packages by `CITYSCAPES_ROOT`. The required directory structure is as follows:
   - `CITYSCAPES_ROOT`
      - `leftImg8bit`
         - `train`
         - `val`
         - `test`
      - `rightImg8bit`
         - `train`
         - `val`
         - `test`
      - `disparity`
         - `train`
         - `val`
         - `test`
      - `camera`
         - `train`
         - `val`
         - `test`  
   The lower levels of the directory structure correspond to city directories and files therein (e.g. `leftImg8bit/test/berlin/berlin_000362_000019_leftImg8bit.png`) and are omitted above for brevity.
3. Create a symbolic link to `CITYSCAPES_ROOT` in the [`data/`](data) directory of the repository and name it `Cityscapes`. In Linux, supposing that `FOG_SIMULATION_ROOT` points to the directory into which you have cloned this repository, this can be performed with
   ```
   cd ${FOG_SIMULATION_ROOT}/data
   ln -s ${CITYSCAPES_ROOT} Cityscapes
   ```
4. Open MATLAB. Change its current folder to the directory [`source/Fog_simulation/Experiments/`](source/Fog_simulation/Experiments) of the cloned repository with
   ```
   cd(FOG_SIMULATION_ROOT);
   cd(fullfile('source', 'Fog_simulation', 'Experiments'));
   ```
   Run the experiment for generating **Foggy Cityscapes-refined** by issuing in the Command Window
   ```
   Cityscapes_trainval_refined_stereogf_beta_0_01_serial;
   ```
   This should create the directory `output/Foggy_Cityscapes/` and populate it with synthetic foggy images in `leftImg8bit_trainval_refined_stereogf_beta_0.01_foggy` as well as corresponding estimated transmittance maps in `leftImg8bit_trainval_refined_stereogf_beta_0.01_transmittance` and depth maps in `depth_stereoscopic_trainval_refined`. **Note**: the results of this experiment will occupy around **6 GB of disk space**.  
   *Foggy Cityscapes-refined* is based on a refined list of 550 Cityscapes images (498 `train` plus 52 `val`) that yield high-quality synthetic foggy images; details are given in our [publication][project_page].

The experiment of step 4 uses a single thread and thus runs for around one day on an Intel Core i7 machine with 16 GB RAM and MATLAB 2017b. However, the implementation of the core MATLAB function [`source/Fog_simulation/Fog_simulation_Cityscapes.m`](source/Fog_simulation/Fog_simulation_Cityscapes.m) of our fog simulation experiments allows **parallel execution** of the experiment. To this end, we provide a few example `bash` scripts in the [experiments directory](source/Fog_simulation/Experiments) which launch experiments on a [Grid Engine](https://en.wikipedia.org/wiki/Oracle_Grid_Engine) cluster for faster execution. These scripts along with the [singled-threaded MATLAB script](source/Fog_simulation/Experiments/Cityscapes_trainval_refined_stereogf_beta_0_01_serial.m) can be consulted for creating a user-specific script for parallel execution depending on the features of the user's system.

Last but not least, to run fog simulation on the `train_extra` split of Cityscapes, the packages `leftImg8bit_trainextra.zip`, `rightImg8bit_trainextra.zip`, `disparity_trainextra.zip` and `camera_trainextra.zip` must be downloaded. Make sure that the extracted directories obey the aforementioned structure, for example
- `CITYSCAPES_ROOT`
  - `leftImg8bit`
     - `train`
     - `val`
     - `test`
     - `train_extra`


### References

1. M. Cordts, M. Omran, S. Ramos, T. Rehfeld, M. Enzweiler, R. Benenson, U. Franke, S. Roth, and B. Schiele: **The Cityscapes Dataset for Semantic Urban Scene Understanding**. In CVPR (2016).
2. L. Wang, H. Jin, R. Yang, and M. Gong: **Stereoscopic Inpainting: Joint Color and Depth Completion from Stereo Images**. In CVPR (2008).
3. R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua, and S. Süsstrunk: **SLIC Superpixels Compared to State-of-the-Art Superpixel Methods**. IEEE Transactions on Pattern Analysis and Machine Intelligence 34(11), 2274-2282 (2012).
4. K. He, J. Sun, and X. Tang: **Guided Image Filtering**. IEEE Transactions on Pattern Analysis and Machine Intelligence 35(6), 1397–1409 (2013).
5. K. He, J. Sun, and X. Tang: **Single Image Haze Removal Using Dark Channel Prior**. IEEE Transactions on Pattern Analysis and Machine Intelligence 33(12), 2341–2353 (2011).
6. K. Tang, J. Yang, and J. Wang: **Investigating Haze-Relevant Features in a Learning Framework for Image Dehazing**. In CVPR (2014).


### Contact

Christos Sakaridis  
csakarid[at]vision.ee.ethz.ch  
http://people.ee.ethz.ch/~csakarid/SFSU_synthetic

[project_page]: <http://people.ee.ethz.ch/~csakarid/SFSU_synthetic/>
[cityscapes]: <https://www.cityscapes-dataset.com/>
[cityscapes_citation]: <https://www.cityscapes-dataset.com/citation/>
[slic_citation]: <https://ivrl.epfl.ch/research/superpixels/>
