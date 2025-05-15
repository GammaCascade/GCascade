# $\gamma$-Cascade
This is a gamma-ray propagation package for Wolfram Mathematica created by Carlos Blanco and Antonio Capanema. $\gamma$-Cascade uses a semi-analytic  approach to model gamma-ray propagation through cosmological distances accounting for attenuation, the formation of electromagnetic cascades,and cosmological redshifting. Our latest paper that explains the inner workings of this package is also the most up to date documentation: https://arxiv.org/abs/2408.03995.

Due to LFS file size restrictions, functional libraries must be downloaded from the [Zenodo repo](https://zenodo.org/doi/10.5281/zenodo.13154969). 

Auxiliary libraries can also be found in the Zenodo repo. Auxiliary libraries are not required to run $\gamma$-Cascade, and are provided for references. See the [documentation](https://arxiv.org/abs/2408.03995) for details. 

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

# V4.1

As of May 15, 2025 $\gamma$-Cascade V4 has been updated to be compatible with Mathematica versions 12.2 and earlier. We have also added to the tutorial notebook a section on how to use $\gamma$-Cascade for obtaining gamma-ray fluxes from cosmological dark matter decay/annihilation.

# V4

As of August 6, 2024 a new version of $\gamma$-Cascade is available, GCascadeV4. V4 implements an assortment of the most widely used EBL models, significantly improves computational precision, and provides new core functionality. Additionally, there is a new method to estimate the uncertainty due to the EBL model. New Libraries have also been created and uploaded to zenodo under the DOI:[10.5281/zenodo.13154970](https://doi.org/10.5281/zenodo.13154970). 
