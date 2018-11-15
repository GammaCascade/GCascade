# GCascade
gamma-ray propagation package for Wolfram Mathematica created by Carlos Blanco

Documentation: https://arxiv.org/abs/1804.00005


# V3
As of November 14, 2018 a new version of GCascade is available, GCascadeV3. Version three drastically improves library loading speeds as well as computational speeds. New Libraries have also been created for this version and consolidated into a handful of .mat files.

New Features:
-Stepsize is not longer required as an input.
-Evolving source spectra now supported for calculation of diffuse gamma-ray background.


# V2
As of June 6, 2018 a new version of GCascade is available, GCascadeV2. Version two drastically improves accuracy at energies below around one PeV. New Libraries have also been created for this version.


 

To use GCascade dowload the .wl and .nb files and place them in the same directory as the libraries folder.

NOTE: Decompressing the 7zip file will produce a single directory titled "LibrariesV3". Inside this directory, there is a set of .mat and .csv files needed to run GCascade. 

The files found in LibrariesV3 should be as follows:

	DiffCascadeSpectralTable.mat
	IMFP.csv
	NormalizedCascadedSecondaryGammaSpectralTable.mat
	NormalizedPairProductionSpectralTable.mat
	testzArray.mat
	zRegIndexArray.mat
