To conserve space, the traces simulated for each transmit are Fourier transformed, and the scripts in this folder (the [visualization](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/tree/main/visualization) folder) can be used to visualize the time-domain data for each transmit. 

The [visualizeTimeDomain.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/visualization/visualizeTimeDomain.m) script convert the Fourier-transformed data for a single transmit simulation back into the time-domain and loops the signal traces through each row of receivers.

The [resamplingTest.m](https://github.com/rehmanali1994/3D-FWI-MultiRowRingArrayUST/blob/main/visualization/resamplingTest.m) script shows the impact that **frequency-domain resampling** can have on the time-domain signal trace.  Basically, if all forward and inverse Fourier transforms are done correctly, the resampling may be done as long as `df` is small enough (so as not to alias the time-domain signal).

c
