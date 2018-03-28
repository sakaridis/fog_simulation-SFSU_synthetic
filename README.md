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
5. Beyond the demo: fog simulation at dataset scale
6. Extra downloads: **CNN models** fine-tuned on *Foggy Cityscapes*


### Requirements

1.  **MATLAB**  
    The code has been developed and tested in MATLAB releases 2016b and 2017b. We therefore recommend using release 2016b or later, noting that these releases are only available for 64-bit systems. If such a configuration is not possible, using an earlier MATLAB release is the recommended (though not tested) alternative.
2.  **C compiler**  
    Users will have to build a binary MEX file for SLIC superpixels themselves in MATLAB (see [instructions](#installation-for-running-the-demo) below), which requires a MATLAB-supported C compiler.


### Installation for running the demo

Our fog simulation pipeline makes use of [SLIC superpixels][slic_citation], for which the algorithm is implemented in the form of a [source MEX file](source/external/SLIC_mex/slicmex.c) for usage in MATLAB. In short, a C compiler has to be configured with MATLAB and then used to build the binary MEX file for SLIC which is called in the pipeline. This [MATLAB doc page](https://www.mathworks.com/help/matlab/matlab_external/what-you-need-to-build-mex-files.html) can serve as a reference.

Steps:
1. Clone this repository with
   ```
   git clone https://github.com/sakaridis/fog_simulation-SFSU_synthetic.git
   ```
2. Make sure a MATLAB-supported C compiler is installed in your system.
3. Open MATLAB and type in the Command Window
   ```
   mex -setup
   ```
   which will guide you through the process of configuring MATLAB's `mex` command with the compiler from step 2. Once this setup has been completed successfully, issuing the above command for a second time should generate a message similar to
   ```
   MEX configured to use 'gcc' for C language compilation.
   ```
   For further details, please consult the [mex command documentation](https://www.mathworks.com/help/matlab/ref/mex.html).
4. Assign the MATLAB variable `fog_simulation_root` with the path to the directory into which you have cloned this repository, by issuing in the Command Window something like
   ```
   fog_simulation_root = '/home/some_user/fog_simulation-SFSU_synthetic';
   ```
5. Change MATLAB's current folder to `fog_simulation_root`, e.g. with
   ```
   cd(fog_simulation_root);
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
cd(fullfile(fog_simulation_root, 'source'));
Demo_fog_simulation_Cityscapes;
```

The demo runs our fog simulation on an example clear-weather image from Cityscapes and writes the results (synthesized foggy image, estimated transmittance map and depth map) under the directory `output/demos/`. The running time for a single image is around 2-3 minutes on an Intel Core i7 machine with 16 GB RAM.


### Foggy Cityscapes-refined

A refined list of 550 Cityscapes images (498 `train` plus 52 `val`) that yield high-quality synthetic foggy images is provided in the file `foggy_trainval_refined_filenames.txt`. Details on the refinement criteria are given in our [publication][project_page].


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
