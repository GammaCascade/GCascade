# $\gamma$-Cascade
This is a gamma-ray propagation package for Wolfram Mathematica created by Carlos Blanco and Antonio Capanema. $\gamma$-Cascade uses a semi-analytic  approach to model gamma-ray propagation through cosmological distances accounting for attenuation, the formation of electromagnetic cascades,and cosmological redshifting.

Due to LFS file size restrictions, functional libraries must be downloaded from the [Zenodo repo](https://zenodo.org/doi/10.5281/zenodo.13154969). 

Auxiliary libraries can also be found in the Zenodo repo. Auxiliary libraries are not required to run $\gamma$-Cascade, and are provided for references. See the documentation for details. 

# Installation

### Step 1

Clone $\gamma$-Cascade repo.

````
git clone https://github.com/GammaCascade/GCascade.git
````

### Step 2
Download LibrariesV4 from the [Zenodo repo](https://zenodo.org/doi/10.5281/zenodo.13154969).

### Step 3
Unzip the LibrariesV4.zip directory into the /GCascade directory. The structure should look as follows:
````
/GCascade
    |-GCascadeV4.wl
    |-LICENSE
    |-README.md
    |-Tutorial.nb
    |-LibrariesV4
        |-...
	|-...
````

### Step 4
Open the Tutorial.nb notebook in Mathematica to learn basic usage.

# V4

As of August 6, 2024 a new version of $\gamma$-Cascade is available, GCascadeV4. V4 implements an assortment of the most widely used EBL models, significantly improves computational precision, and provides new core functionality. Additionally, there is a new method to estimate the uncertainty due to the EBL model. New Libraries have also been created and uploaded to zenodo under the DOI:[10.5281/zenodo.13154970](https://doi.org/10.5281/zenodo.13154970). 


# V3

Originl paper and documentation for versions up to V3: https://arxiv.org/abs/1804.00005

As of November 14, 2018 a new version of GCascade is available, GCascadeV3. V3 drastically improves library loading speeds as well as computational speeds. New Libraries have also been created for this version and consolidated into a handful of .mat files.

NOTE: Since LibrariesV3 is held in git LFS, you might need to download the raw archive directly if it isn't downloaded correctly during cloning. (git clone sometimes downloads a file named "LibrariesV3" without any contents) 

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
