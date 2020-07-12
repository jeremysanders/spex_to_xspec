# spex_to_xspec

spex_to_xspec is a Python script to generate an APEC-format table model (http://www.atomdb.org/) from the SPEX spectral modelling and fitting code  (https://www.sron.nl/astrophysics-spex). This allows the user to load a SPEX-generated plasma model into Xspec (https://heasarc.gsfc.nasa.gov/lheasoft/xanadu/xspec) to use the same environment for both APEC and SPEX spectral fitting.

## Getting Started

Download the spex_to_xspec.py script. You will need to setup the SPEX environment.

There are some adjustable parameters you can modify at the top of the script.  These include the default filename root part (to be used in Xspec when loading the APEC table), the list of temperatures to tabulate the SPEX model for and the continuum energy binning.

When you are ready, run the script using python.

By default the temperature grid is the same as the APEC model. There is a higher resolution temperature grid (0.01 to 100 keV in 201 log bins) option which is slower. The output files are currently very large (1GB+), as every line in the energy range (by default 0.05 to 15 keV in 300 log bins) is dumped. You can compress the output files with gzip and Xspec is able to read them.

Note that the program dumps out a lot of data into temporary files (called tmp_*) which are later interpreted. If you stop the program mid way through, you will need to clean up these files.

To use the table in xspec use
```
xset APECROOT spex
```
Where `spex` is the root specified at the top of the script.

### Prerequisites

 1. Python 2.7 or Python 3.3+
 2. Astropy (http://www.astropy.org/)
 3. SPEX 3.04 (https://www.sron.nl/astrophysics-spex). Future versions may work.

## How it works

The script generates a line list and a continuum for each of the temperature grid points in ASCII format from SPEX. It then interprets these and converts them to FITS format. Generating the line list file is fairly simple, but to calculate the continuum for each element, the continuaa are dumped with very high metallicity so that the hydrogen and helium continuaa can be subtracted with good numerical precision.

By default the script switches to the latest SPEX plasma model (`var calc new`).

## Future improvements

 1. Store very weak lines in the pseudo-continuum to reduce the file size
 2. Do not store the temporary files in the current directory
 3. Clean up temporary files on incomplete run

## Authors

* **Jeremy Sanders**  - [jeremysanders](https://github.com/jeremysanders)

## License

This project is licensed under the MIT License - see the [LICENSE.txt](LICENSE.txt) file for details
