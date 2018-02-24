# OCT-B-Scan-Analysis
This MATLAB program analyzes a 2-configuration OCT B-Scan experiment with multiple participants. The program assumes image output from Heidelberg Engineering software and requires each image to be given both in segmented and unsegmented form.

All MATLAB files are given both in .m format and in printed HTML format. Further documentation is given within the code.

## Usage
To perform an automatic analysis of an entire experiment, run the script processExperiment. The following should be modified:
* analysisDir - change to desired directory for output
* octNoSeg - change to desired directory for unsegmented versions of the input
* octSegDir - change to desired directory for segmented versions of the input
* numOfParticipants - change to number of participants in the experiment

The program assumes that within each input directory, there are numOfParticipants directories labeled as 1,2...,numOfParticipants. Within each of the latter, there should be two directories, one for each configuration labled in '####i' format, where i=1,2. Within each of the latter, there should be a constant number of OCT B-Scan images.

By default, the program removes retinal blood vessel signals from the image prior to computation of statistics. This could be altered by modifying the function computeStats in evaluateB.m, as described within the code.

## Reference
This program makes use of the functions findpeaksG and finsqaurepulse by Thomas C. O'Haver (2014). See [Interactive Signal Processing Tools](https://terpconnect.umd.edu/~toh/spectrum/SignalProcessingTools.html) for further information on these functions.

## Citing
This program is part of "Assessing the Effect of Contact Lenses on the Image Quality of Retinal Spectral Domain Optical Coherence Tomography Using Automated Image Analysis" by Yinon Shapira, Talia Aviram, Omer Granak, Igor Viner, Erez Ribak, Eitan Z Blumenthal (2018). The manuscript has been submitted to peer-review.
