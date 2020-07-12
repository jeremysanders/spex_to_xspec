#!/usr/bin/env python3

# Program to dump out spex mekal model at various temperature grid points
# The lines and continuum are then gathered up, and stored in an apec-format
# table model for xspec

# Jeremy Sanders 2005-2018

# Version 2.0 (2018-07-31): Support SPEX 3.04

import os
import subprocess
import os.path

import numpy as np
from astropy.io import fits

#####################################################################
# Adjustable parameters

# output root (using spex version) for apec format filenames
# creates outroot_(line|coco).fits
outroot = 'spex'

##########
# temperature grid: uncomment line and comment other to select

# default APEC temperature grid
#temperatures = np.logspace(np.log10(0.0008617385), np.log10(86.17385),51)

# increased numbers of sample points between 0.01 and 100 keV
temperatures = np.logspace(np.log10(0.01), np.log10(100), 201)

# for testing
#temperatures = np.array([1,2])
##########

# energy range and stepping to sample continuum (log spacing used)
contminenergy = 0.05
contmaxenergy = 15.
contenergysteps = 300

# energy range and stepping to sample pseudo-continuum (log spacing used)
pcontminenergy = 0.05
pcontmaxenergy = 15.
pcontenergysteps = 2048

# lines lower than this flux (photon cm^3/s) are put into a
# pseudo-continuum rather than stored separately.

# The APEC default is 1e-20, but this produces many fewer lines using
# this for SPEX
minepsilon = 1e-22

# where to put output files
tmpdir = 'workdir'

###########################################################

# location of spex installation
try:
    spexroot = os.environ['SPEX90']
except KeyError:
    raise RuntimeError('Please set SPEX90 and initialize the SPEX environment')

# executable to use to run spex
spexexecutable = os.path.join(spexroot, 'bin/spex')

# for checking for numerical data
digits = '0123456789'

# conversions
keV_K = 11.6048e6
keV_erg = 1.6022e-9

# convert from unit norm to cm3 in spex
norm_factor_cm3 = 1e58

# identify element names with element numbers in spex
elements = (
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O',
    'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
    'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti',
    'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn')

# create dict to convert element names into numbers
element_nums = {}
for num, element in enumerate(elements):
    element_nums[element] = num+1

# these are the apec elements
apec_elements = [
    'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Al', 'Si', 'S', 'Ar',
    'Ca', 'Fe', 'Ni']

# with hydrogen
all_elements = ['H']+apec_elements

# roman numerals (bah)
roman_numerals = (
    'I', 'II', 'III', 'IV', 'V', 'VI', 'VII',
    'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV',
    'XV', 'XVI', 'XVII', 'XVIII', 'XIX', 'XX',
    'XXI', 'XXII', 'XXIII', 'XXIV', 'XXV', 'XXVI', 'XXVII',
    'XXVIII', 'XXIX', 'XXX')

# dict to convert numerals to numbers
roman_to_number = {}
for num, numeral in enumerate(roman_numerals):
    roman_to_number[numeral] = num+1

# abundance to use to get continuua of elements more exactly
continuum_mult = 1000.

class Line:
    def __init__(self, element, ion, wavelength, epsilon, energy):
        self.element = element
        self.ion = ion
        self.wavelength = wavelength
        self.epsilon = epsilon
        self.energy = energy

def deleteFile(f):
    """For debugging."""
    #os.unlink(f)

def writeScriptElements(fobj, elements, val):
    """Write commands to set abundance to val."""
    for e in elements:
        print('par %02i val %f' % (element_nums[e], val), file=fobj)

def writeScriptLines(fobj, T):
    """Write script to output lines."""
    # set temperature
    print('par t val %e' % T, file=fobj)

    # switch on all apec elements
    writeScriptElements(fobj, apec_elements, 1)

    # compute model
    print('calc', file=fobj)

    # dump out lines
    outfile = os.path.join(tmpdir, 'tmp_lines_T%010f' % T)
    print('ascdump file %s 1 1 line' % outfile, file=fobj)

    writeScriptElements(fobj, apec_elements, 0)

def writeScriptContinuua(fobj, T):
    # set temperature
    print('par t val %e' % T, file=fobj)

    print('calc', file=fobj)
    outfile = os.path.join(tmpdir, 'tmp_conti_T%010f_%s' % (T, 'H'))
    print('ascdump file %s 1 1 tcl' % outfile, file=fobj)

    for el in apec_elements:
        writeScriptElements(fobj, (el,), continuum_mult)
        print('calc', file=fobj)
        outfile = os.path.join(tmpdir, 'tmp_conti_T%010f_%s' % (T, el))
        print('ascdump file %s 1 1 tcl' % outfile, file=fobj)
        writeScriptElements(fobj, (el,), 0)

def writeScript(fobj, T):
    """Write header to script sent to spex."""

    print('egrid log %e:%e %i' % (
        contminenergy, contmaxenergy,
        contenergysteps), file=fobj)
    # switch to the latest spex model
    print('var calc new', file=fobj)
    print('var newmekal all true', file=fobj)
    print('comp cie', file=fobj)
    print('abundance ag', file=fobj)

    writeScriptElements(fobj, elements[1:], 0)
    writeScriptLines(fobj, T)
    writeScriptContinuua(fobj, T)
    print('quit', file=fobj)

def generateOutput():
    """Process all the temperatures and generate output files."""

    for T in temperatures:
        fname = os.path.join(tmpdir, 'tmp_spex_T%010f.script' % T)
        with open(fname, 'w') as fout:
            writeScript(fout, T)

        with open(fname) as fin:
            subprocess.call([spexexecutable], stdin=fin)
        deleteFile(fname)

def makeLineHDU(lines, T, totflux):
    """Given lines list, produce line HDU.

    lines is (element, ion, wavelength, epsilon, energy) list."""

    # sort lines by element and ion and energy
    lines.sort(key=lambda x: (x.element, x.ion, 1/x.wavelength))

    # construct up FITS table to APEC format
    col_lambda = fits.Column(
        name='Lambda', format='1E', unit='A',
        array=[i.wavelength for i in lines])
    col_lambda_err = fits.Column(
        name='Lambda_Err', format='1E', unit='A',
        array=np.zeros( (len(lines),) ) + np.nan)
    col_epsilon = fits.Column(
        name='Epsilon', format='1E',
        unit='photons cm^3 s^-1',
        array=[v.epsilon for v in lines])
    col_epsilon_err = fits.Column(
        name='Epsilon_Err', format='1E',
        unit='photons cm^3 s^-1',
        array=np.zeros((len(lines),)) + np.nan)
    col_element = fits.Column(
        name='Element', format='1J',
        array=[i.element for i in lines])
    col_ion = fits.Column(
        name='Ion', format='1J',
        array=[i.ion for i in lines])

    col_upperlev = fits.Column(
        name='UpperLev', format='1J',
        array=np.zeros((len(lines),)) + 2)
    col_lowerlev = fits.Column(
         name='LowerLev', format='1J',
         array=np.zeros((len(lines),)) + 1)

    tabhdu = fits.BinTableHDU.from_columns([
        col_lambda, col_lambda_err, col_epsilon,
        col_epsilon_err, col_element, col_ion,
        col_upperlev, col_lowerlev])

    tabhdu.name = 'EMISSIVITY'
    h = tabhdu.header
    h['HIERARCH TEMPERATURE'] = T*keV_K
    h['XTEMP'] = T

    # fixme wrong below (erg not photon)
    h['TOT_LINE'] = totflux
    h['N_LINES'] = len(lines)
    return tabhdu

def interpretDumpedLines(T):
    """Interpret dumped lines file.

    Returns new HDU, number of lines, and number of elements
    """

    print('Interpreting dumped lines for T=%g' % T)

    totflux = 0.
    elements = set()
    weak_lines = []
    lines = []
    outfile = os.path.join(tmpdir, 'tmp_lines_T%010f.asc' % T)
    for line in open(outfile):
        p = line.strip().split()
        # numerical line
        if p[0][0] in digits:
            # horribly, have to hard code in column numbers here
            element = element_nums[line[9:12].strip()]
            ion = roman_to_number[ line[12:17].strip() ]
            wavelength = float( line[102:115] )
            energy_keV = float( line[87:100] )

            # convert from total photon flux to normalised photon flux
            epsilon = float( line[117:126] ) / norm_factor_cm3

            # skip lines out of energy range
            if energy_keV<contminenergy or energy_keV>contmaxenergy:
                continue

            line = Line(element, ion, wavelength, epsilon, energy_keV)

            if epsilon > minepsilon:
                # keep track of total flux in ergs
                totflux += energy_keV*keV_erg*epsilon

                elements.add(element)
                lines.append(line)
            else:
                weak_lines.append(line)

    deleteFile(outfile)

    print('T=%g, %i strong lines, %i weak lines' % (
        T, len(lines), len(weak_lines)))

    tabhdu = makeLineHDU(lines, T, totflux)
    return tabhdu, len(lines), len(elements), weak_lines

def readContinuum(filename):
    """Take spex dumped model spectrum file, and extract continuum."""
    outenergy = []
    outval = []
    for line in open(filename):
        p = line.strip().split()
        if p[0][0] in digits:
            outenergy.append(float(p[1]))
            outval.append(float(p[2]) * 1e44 / norm_factor_cm3)
    return (np.array(outenergy), np.array(outval))

def interpretDumpedContinuum(T):
    """Interpret dumped continuum file.
    """

    print('Interpreting dumped continuum for T=%g' % T)

    # read in continum from each file
    continuua = {}

    for element in all_elements:
        filename = os.path.join(
            tmpdir, 'tmp_conti_T%010f_%s.asc' % (T, element))
        energy, vals = readContinuum(filename)
        continuua[element] = vals
        deleteFile(filename)

    # subtract H continuum from each element except hydrogen
    # also divide by abundance these were generated at
    for element in apec_elements:
        continuua[element] -= continuua['H']
        continuua[element] /= continuum_mult

    # construct table
    contformat = '%iE' % contenergysteps
    col_element = fits.Column(
        name='Z', format='1J',
        array=[element_nums[i] for i in all_elements])
    col_rmJ = fits.Column(
        name='rmJ', format='1J', array=np.zeros(len(all_elements)))
    col_N_Cont = fits.Column(
        name='N_Cont', format='1J',
        array= [contenergysteps]*len(all_elements))
    col_E_Cont = fits.Column(
        name='E_Cont', format=contformat,
        unit='keV', array=np.resize(
            energy, (len(all_elements), contenergysteps)))

    col_Continuum = fits.Column(
        name='Continuum', format=contformat,
        unit='photons cm^3 s^-1 keV^-1',
        array=[continuua[i] for i in all_elements])
    col_Cont_Err = fits.Column(
        name='Cont_Err', format=contformat,
        unit='photons cm^3 s^-1 keV^-1',
        array=np.zeros( (len(all_elements), contenergysteps) ))

    # create zero pseudo-continuum to fill in later
    pcontformat = '%iE' % pcontenergysteps
    col_N_Pseudo = fits.Column(
        name='N_Pseudo', format='1J',
        array= [pcontenergysteps]*len(all_elements))
    col_E_Pseudo = fits.Column(
        name='E_Pseudo', format=pcontformat,
        unit='keV', array=np.resize(
            energy, (len(all_elements), pcontenergysteps)))
    col_Pseudo = fits.Column(
        name='Pseudo', format=pcontformat,
        array=np.zeros( (len(all_elements), pcontenergysteps) ))
    col_Pseudo_Err = fits.Column(
        name='Pseudo_Err', format=pcontformat,
        array=np.zeros( (len(all_elements), pcontenergysteps) ))

    tabhdu = fits.BinTableHDU.from_columns([
        col_element, col_rmJ,
        col_N_Cont, col_E_Cont,
        col_Continuum, col_Cont_Err,
        col_N_Pseudo, col_E_Pseudo,
        col_Pseudo, col_Pseudo_Err])

    tabhdu.name = 'EMISSIVITY'
    h = tabhdu.header
    h['HIERARCH TEMPERATURE'] = T*keV_K
    h['XTEMP'] = T
    h['DENSITY'] = 1.0
    # sum flux
    totcoco = 0.
    for i in all_elements:
        totcoco += continuua[i].sum()
    h['TOT_COCO'] = totcoco

    return (
        tabhdu,
        len(all_elements),
        contenergysteps*len(all_elements),
        pcontenergysteps*len(all_elements)
    )

def interpretAllLines():
    """Interpret spex dumped spectra."""

    # generate HDUs for each temperature
    hdus = []
    Nelement = []
    Nline = []
    weaklinelist = []
    for T in temperatures:
        hdu, numlines, numelements, weaklines = interpretDumpedLines(T)
        Nline.append(numlines)
        Nelement.append(numelements)
        hdus.append(hdu)
        weaklinelist.append(weaklines)

    # construct HDU describing parameters
    col_kT = fits.Column(
        name='kT', format='1E', unit='keV',
        array=temperatures)
    col_EDensity = fits.Column(
        name='EDensity', format='1E', unit='cm**-3',
        array=np.zeros(len(temperatures))+1)
    col_Nelement = fits.Column(
        name='Nelement', format='1J',
        array=Nelement)
    col_Nline = fits.Column(
        name='Nline', format='1J',
        array=Nline)
    tabhdu = fits.BinTableHDU.from_columns([
        col_kT, col_EDensity, col_Nelement, col_Nline])
    tabhdu.name = 'PARAMETERS'

    # make output file containing all lines
    hdulist = fits.HDUList([fits.PrimaryHDU(), tabhdu] + hdus)
    hdulist.writeto('%s_line.fits' % outroot, overwrite=True)

    return weaklinelist

def computePseudoContinuum(hdu, weaklines):
    """Compute pseudo continuum and enter into continuum HDU."""

    energyedges = np.logspace(
        np.log10(pcontminenergy), np.log10(pcontmaxenergy),
        pcontenergysteps+1)

    for i, element in enumerate(all_elements):
        elidx = element_nums[element]
        energies = [line.energy for line in weaklines if line.element==elidx]
        epsilons = [line.epsilon for line in weaklines if line.element==elidx]

        summedlines, edgesout = np.histogram(
            energies, weights=epsilons, bins=energyedges)
        # divide by bin width to convert to photon cm^3/s/keV
        flux = summedlines / (energyedges[1:]-energyedges[:-1])
        hdu.data.field('Pseudo')[i,:] = flux

def interpretAllContinuum(weaklinelist):
    """Build up continuum output file."""

    # make continuum HDUs for each temperature
    hdus = []
    NElement = []
    NCont = []
    NPseudo = []
    for T, weaklines in zip(temperatures, weaklinelist):
        hdu, numelem, numcont, numpseudo = interpretDumpedContinuum(T)
        computePseudoContinuum(hdu, weaklines)
        hdus.append(hdu)
        NElement.append(numelem)
        NCont.append(numcont)
        NPseudo.append(numpseudo)

    # construct HDU describing parameters
    col_kT = fits.Column(
        name='kT', format='1E', unit='keV',
        array=temperatures)
    col_EDensity = fits.Column(
        name='EDensity', format='1E', unit='cm**-3',
        array=np.zeros(len(temperatures))+1)
    col_NElement = fits.Column(
        name='NElement', format='1J',
        array=NElement)
    col_NCont = fits.Column(
        name='NCont', format='1J',
        array=NCont)
    col_NPseudo = fits.Column(
        name='NPseudo', format='1J',
        array=NPseudo)
    tabhdu = fits.BinTableHDU.from_columns([
        col_kT, col_EDensity, col_NElement, col_NCont, col_NPseudo
    ])
    tabhdu.name = 'PARAMETERS'

    # make output file containing the continuum
    hdulist = fits.HDUList([fits.PrimaryHDU(), tabhdu] + hdus)
    hdulist.writeto('%s_coco.fits' % outroot, overwrite=True)

def main():
    """Main routine."""
    generateOutput()
    weaklinelist = interpretAllLines()
    interpretAllContinuum(weaklinelist)

if __name__ == '__main__':
    main()
