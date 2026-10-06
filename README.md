# S-bursts from T8 brown dwarf WISE J0623-0456

## Description

Custom codes that went into creating plots for article (under review) "A cold brown dwarf with a dense magnetosphere". They include:

[1] Example scripts for visualization of the dynamic spectra and S-burst property histograms

[2] Example scripts for calculating the Hill-Pontius radius and the Alfvén radius

The scripts offer a quick approach to reproduce the results in the article, for instance plotting the dynamic spectra of the S-bursts. Processed dynamic spectra and the S-burst properties could be accessed at https://disk.pku.edu.cn/link/AA9ABB9B2D79C249B39DAB270C29DBD33C.

The complete pipeline for reducing FAST observation data have been published in "Starspots as the origin of ultrafast drifting radio bursts from an active M dwarf" (https://www.science.org/doi/10.1126/sciadv.adw6116) and the associated Zenodo (https://doi.org/10.5281/zenodo.15352945) and Github links (https://github.com/jiale-radioastro/ADLeo_ultrafast_drifting_bursts). The raw FAST observational data used in this work are available from the FAST archive (http://fast.bao.ac.cn, project ID: PT2024_0017, PT2025_0037) or reasonable requests to the authors. Note the total amount of raw data of these two projects is around 30 TB.

## Requirements

The scripts were run using Python 3.9.23 with the following package versions:

- NumPy 1.26.4
- Matplotlib 3.5.2
- Astropy 6.0.1
- SciPy 1.13.1

Other Python and package versions are likely to work as well, but have not been tested. 
All dependencies can be installed using `pip` or `conda`. Please refer to each package's documentation for installation instructions.


