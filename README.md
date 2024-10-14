# 3D-FWI-MultiRowRingArrayUST
3D Full-Waveform Inversion for a Multi-Row Ring-Array UST System SImulated in k-Wave

Ultrasound tomography (UST) is a medical imaging system that uses the transmission of ultrasound through tissue to create images of the speed of sound.  Currently, UST is based on a single-row elevation-focused ring-array transducer.  This ring-array is translated vertically to reconstruct a series of 2D slices that are assembled into a 3D volume.  The main algorithm used to reconstruct these 2D slices is 2D frequency-domain full-waveform inversion (FWI): see our previous work on 2D slice-wise FWI ([rehmanali1994/WaveformInversionUST](https://github.com/rehmanali1994/WaveformInversionUST)).  

In this work, we extend the ring-array geometry to multiple rows of elements.  By emitting cylindrical wave transmits from this multi-row ring array, we can keep the number of transmisssions to a minimum while still capturing the full 3D insonification needed to reconstruct the full volume.  With this imaging geometry, 2D FWI can be applied to each row of receive elements to reconstruct an image slice for each row.  Here, we compare 3D FWI to 2D slicewise FWI in numerical breast phantoms.  

We provide sample data and algorithms presented in

> Ali, R., Jin, G., Singh, M., Mitcham, T., & Duric, N. (2024, _In Review_). 3D Frequency-Domain Full Waveform Inversion for Whole-Breast Imaging with a Multi-Row Ring Array. _IEEE Open Journal of Ultrasonics, Ferroelectrics, and Frequency Control._

You can also reference a static version of this code by its DOI number:
[![DOI](https://zenodo.org/badge/871745761.svg)](https://doi.org/10.5281/zenodo.13924218)

