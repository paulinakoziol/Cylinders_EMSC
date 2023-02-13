# Cylinders_EMSC

### Resonant Mie-type Scattering Extended Multiplicative Signal Correction for samples with cylindrical domains ###

Author: Paulina Koziol

Development of this code was supported by the National Science Centre, Poland (“Improving 3D macromolecule orientation determination based on polarized IR chemical imaging by optimization of scattering removal algorithms”, Grant No. 2019/35/N/ST4/02481).

## Background ##
 
This GPU-based algorithm is used to correct Mie-type scattering artifacts observed during IR transmission measurements of samples with cylindrical domains. Due to the cylindrical symmetry of the samples, a case o linearly polarized light is also considered apart from an unpolarized light case. An exemplary data of PCL fiber are provided with the code. This algorithm also requires a reference spectrum that is similar to the sample under investigation. PCL reference spectrum is provided.  

A research paper describing this algorithm and presenting examples of its usage is currently being published. 

The fundamental physics is described in [1]. 

The algorithm is based on the RMieS-EMSC algorithm described in [2] and [3], with an open source code available: https://github.com/GardnerLabUoM/RMieS 

[1] I. L. Rasskazov, R. Singh, P. S. Carney, and R. Bhargava, “Extended Multiplicative Signal Correction for Infrared Microspectroscopy of Heterogeneous Samples with Cylindrical Domains,” Appl. Spectrosc., vol. 73, no. 8, pp. 859–869, 2019, doi: 10.1177/0003702819844528

[2] P. Bassan et al., “RMieS-EMSC correction for infrared spectra of biological cells: extension using full Mie theory and GPU computing.,” J. Biophotonics, vol. 3, no. 8–9, pp. 609–20, Aug. 2010, doi: 10.1002/jbio.201000036

[3] P. Bassan et al., “Resonant Mie scattering (RMieS) correction of infrared spectra from highly scattering biological samples.,” Analyst, vol. 135, no. 2, pp. 268–77, Feb. 2010, doi: 10.1039/b921056c

## License ##
 
This code is licensed under the GNU Lesser General Public License v3.0.

## Installation ##

Download the package, unzip it and open main.m script in MATLAB. This script is ready to be used with model PCL data provided in the package. However, for the acceleration purpose, functions use GPU implementation and MATLAB compatible GPU unit is necessary. 

Scattering free PCL spectrum is used as a build-in reference spectrum. 

Tips to use Cylinders_EMSC package:

1. For the use of the cylinders_EMSC function, prepare your data to be in the following form: 
	- Column vector WN (size: [K x 1]) for the wavenumber values in ascending order. 
	- An IR data matrix called dataR (size: [N x K]) with one spectrum per row. N is the number of spectra to be corrected. 
An example of the proper shaping of the 3D data is provided in the main.m script. 

2. Define correction options. An example for the model PCL fiber correction (unpolarized light case) is provided below: 
    
```matlab
correction_options = [ ...
    898  ;      % 1. Lower wavenumber range
    3845 ;      % 2. Upper wavenumber range
    10   ;      % 3. Number of iterations
    7    ;      % 4. Number of principal components used (7 default)
    5    ;      % 5. Lower range for scattering particle radius / um
    15   ;     % 6. Upper range for scattering particle radius / um
    1.3  ;      % 7. Lower range for average refractive index
    1.7  ;      % 8. Upper range for average refractive index
    10   ;      % 9. Number of values for each scattering parameter (a,b,d) default 10
    1    ;      % 10. Gram-Schmidt Process option, 0 = no, 1 = yes. (1 recommended)
    2    ;      % 11. Which reference spectrum, 1 = Matrigel, 2 = PCL spectrum , 3 = Mean spectrum
    pi/2 ;      % 12. Incidence angle (zeta) (perpendicular - 90deg = pi/2)
    -1  ];      % 13. Cylinder orientation with respect to linear polarization (degrees):
                %     from 0 for paralell to 90 for perpendicular orientation,
                %     -1 for unpolarized light.                   
```            

Adjust your correction options as necessary. 

3. cylinders_EMSC function gives the following output: 
	- WN_corr - a wavenumber vector for corrected spectra.
	- correctedSpectra - a matrix with corrected spectra. 
	- history - 3D matrix with history for the corrected spectra. 
	- timing - calculation time for each spectrum. 

4. cylinders_EMSC function has a PCL reference and Matrigel reference (for cells and tissues) spectra implemented depending on the choice in correction_options. Another option is to use the mean spectrum of your data. Also, the user might provide its own reference spectrum using the following command: 
	[WN_out, ZCorr, history, timing] = cylinders_EMSC(WN, ZRaw, correction_options, WN_Ref_in, Ref_in); 
where WN_Ref_in and Ref_in are the wavenumber vector and reference spectrum, respectively.   


Other remarks:

1. The CO<sub>2</sub> region of the spectrum is automatically down-weighted and ignored in this algorithm.

2. The corrected spectra are normalized using the BASIC un-weighted EMSC approach.

