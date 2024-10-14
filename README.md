# 3D-FWI-MultiRowRingArrayUST
3D Full-Waveform Inversion for a Multi-Row Ring-Array UST System SImulated in k-Wave

Ultrasound tomography (UST) is a medical imaging system that uses the transmission of ultrasound through tissue to create images of the speed of sound.  Currently, UST is based on a single-row elevation-focused ring-array transducer.  This ring-array is translated vertically to reconstruct a series of 2D slices that are assembled into a 3D volume.  The main algorithm used to reconstruct these 2D slices is 2D frequency-domain full-waveform inversion (FWI): see our previous work on 2D slice-wise FWI ([rehmanali1994/WaveformInversionUST](https://github.com/rehmanali1994/WaveformInversionUST)).  

In this work, we extend the ring-array geometry to multiple rows of elements.  We simulate a 22-cm diameter ring array with 32 rows of 256 elements placed circumferentially around the ring (2.4 mm between each row).  By emitting cylindrical-wave transmits from this multi-row ring array, we can keep the number of transmisssions to a minimum while still capturing the full 3D insonification needed to reconstruct the full volume.  With this imaging geometry, 2D FWI can be applied to each row of receive elements to reconstruct an image slice for each row.  Here, we compare 3D FWI to 2D slicewise FWI in numerical breast phantoms.  

We provide sample data and algorithms presented in

> Ali, R., Jin, G., Singh, M., Mitcham, T., & Duric, N. (2024, _In Review_). 3D Frequency-Domain Full Waveform Inversion for Whole-Breast Imaging with a Multi-Row Ring Array. _IEEE Open Journal of Ultrasonics, Ferroelectrics, and Frequency Control._

You can also reference a static version of this code by its DOI number:
[![DOI](https://zenodo.org/badge/871745761.svg)](https://doi.org/10.5281/zenodo.13924218)

## k-Wave Simulation

The code used to run the [k-Wave](http://www.k-wave.org/) simulation involves a 2-step process:

1) [GenKWaveSimInfoFromMRI.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/GenKWaveSimInfoFromMRI.m) creates a MAT file that is stored in the [sim_info](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/tree/main/sim_info) folder. This MAT file (either _LeftBreastMRI.mat_ or _RightBreastMRI.mat_) contains all the information (medium, multi-row ring-array geometry, and pulse excitation) needed to simulate the UST system. [GenKWaveSimInfoFromMRI.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/GenKWaveSimInfoFromMRI.m) uses [sampled_cylinder.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/phantoms/sampled_cylinder.m) to create the multi-row ring-array transducer and [BreastPhantomVolumetricMRI.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/phantoms/BreastPhantomVolumetricMRI.m) to produce the material properties used in the simulation. [BreastPhantomVolumetricMRI.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/phantoms/BreastPhantomVolumetricMRI.m) relies on the BreastPhantomFromMRI.mat found under the [releases](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/releases) tab.  Please place BreastPhantomFromMRI.mat into the [phantoms](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/tree/main/phantoms) folder.  Note that [GenKWaveSimInfoFromMRI.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/GenKWaveSimInfoFromMRI.m) can run two different simulations: one corresponding to the left breast (`leftBreast = true`) and the other corresponding to the right breast (`leftBreast = false`).  This option is passed into [BreastPhantomVolumetricMRI.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/phantoms/BreastPhantomVolumetricMRI.m) function where it selects one of the breast for the phantom.
2) After generating the MAT file in [Simulations/sim_info](https://github.com/rehmanali1994/FrequencyDifferencing/tree/main/Simulations/sim_info), we run the actual [k-Wave](http://www.k-wave.org/) simulation using [kWaveSimLauncher.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/kWaveSimLauncher.m).  [kWaveSimLauncher.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/kWaveSimLauncher.m) loops through each single-element transmit (function calls to [GenRFDataSingleTxKWave.m](https://github.com/rehmanali1994/FrequencyDifferencing/blob/main/Simulations/GenRFDataSingleTxKWave.m)). 
 The simulated data for each transmit is then stored in MAT files in the [Simulations/sim_data](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/tree/main/sim_data) folder.  To conserve space, the traces simulated for each transmit are Fourier transformed, and the scripts in the [visualization](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/tree/main/visualization) folder can be used to visualize the time-domain data for each transmit.  At the end of [kWaveSimLauncher.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/kWaveSimLauncher.m), the function call to 
[AssembleSimData.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/AssembleSimData.m) assembles the simulated data from each indivdual transmit/MAT-file in the [Simulations/sim_data](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/tree/main/sim_data) folder into a single MAT file (either _LeftBreastMRI.mat_ or _RightBreastMRI.mat_) containing the 3D UST dataset needed for FWI in the [datasets](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/tree/main/datasets) folder.  Note that the complete time-series data is too large to store in a single file for the full 3D dataset.  Instead, we extract only the frequencies listed in [kWaveSimLauncher.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/kWaveSimLauncher.m) to keep the file size for the full 3D dataset small. 

## Reconstruction Code

Two different scripts are used to run FWI (either 2D slicewise or 3D FWI) and their results are compared:
1) [MultiFrequencyWaveformInvFromData.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/MultiFrequencyWaveformInvFromData.m) is used to run **2D slicewise FWI** and outputs the result into [results/results2D](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/tree/main/results/results2D).
2) [MultiFrequencyWaveformInvFromData3D.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/MultiFrequencyWaveformInvFromData3D.m) is used to run **3D FWI** and outputs the result into [results/results3D](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/tree/main/results/results3D).

The key functions used in the FWI scripts above are: 
1) [solveHelmholtzBornSeries.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/solvers4FWI/solveHelmholtzBornSeries.m) uses the convergent Born series method to

Please cite the following paper for the convergent Born series method:
```BibTeX
@article{saha2018machine,
  title={A machine learning approach to radiogenomics of breast cancer: a study of 922 subjects and 529 DCE-MRI features},
  author={Saha, Ashirbani and Harowicz, Michael R and Grimm, Lars J and Kim, Connie E and Ghate, Sujata V and Walsh, Ruth and Mazurowski, Maciej A},
  journal={British journal of cancer},
  volume={119},
  number={4},
  pages={508--516},
  year={2018},
  publisher={Nature Publishing Group UK London}
}```

