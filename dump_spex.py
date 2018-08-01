#!/usr/bin/env python

# Program to dump out spex mekal model at various temperature grid points
# The lines and continuum are then gathered up, and stored in an apec-format
# table model for xspec

# Jeremy Sanders 2005-2018

# Version 2.0 (2018-07-31): Support SPEX 3.04

from __future__ import print_function, division

import os
import subprocess

import numpy as np
from astropy.io import fits

#####################################################################
# Adjustable parameters

# output root (using spex version) for apec format filenames
# creates outroot_(line|coco).fits
outroot = 'spex'

####
# temperature grid: uncomment line and comment other to select

# default APEC temperature grid
# note: SPEX seems to hang indefinitely with these very low temperatures
#temperatures = np.logspace(np.log10(0.0008617385), np.log10(86.17385),51)

# increased numbers of sample points between 0.01 and 100 keV
temperatures = np.logspace(np.log10(0.01), np.log10(100),201)

# for testing
#temperatures = np.linspace(1,5)
####

# energy range and stepping to sample continuua (log spacing used)
contminenergy = 0.05
contmaxenergy = 15.
contenergysteps = 300

###########################################################

# location of spex installation
try:
    spexroot = os.environ['SPEX90']
except KeyError:
    raise RuntimeError('Set the SPEX90 to point to your SPEX installation')

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
continuum_mult = 10000.

def deleteFile(f):
    """For debugging."""
    os.unlink(f)

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
    outfile = 'tmp_lines_T%010f' % T
    print('ascdump file %s 1 1 line' % outfile, file=fobj)

    writeScriptElements(fobj, apec_elements, 0)

def writeScriptContinuua(fobj, T):
    # set temperature
    print('par t val %e' % T, file=fobj)

    print('calc', file=fobj)
    outfile = 'tmp_conti_T%010f_%s' % (T, 'H')
    print('ascdump file %s 1 1 tcl' % outfile, file=fobj)

    for el in apec_elements:
        writeScriptElements(fobj, (el,), continuum_mult)
        print('calc', file=fobj)
        outfile = 'tmp_conti_T%010f_%s' % (T, el)
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
        fname = 'tmp_spex_T%010f.script' % T
        with open(fname, 'w') as fout:
            writeScript(fout, T)

        with open(fname) as fin:
            subprocess.call([spexexecutable], stdin=fin)
        deleteFile(fname)

def interpretDumpedLines(T):
    """Interpret dumped lines file.

    Returns new HDU, number of lines, and number of elements
    """

    print('Interpreting dumped lines for T=%g' % T)

    totflux = 0.
    elements = {}
    lines = []
    outfile = 'tmp_lines_T%010f.asc' % T
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
            strength = float( line[117:126] ) / norm_factor_cm3

            # keep track of total flux in ergs
            totflux += energy_keV*keV_erg*strength

            # skip lines out of energy range
            if energy_keV<contminenergy or energy_keV>contmaxenergy:
                continue

            lines.append( (element, ion, wavelength, strength) )
            elements[element] = True
    deleteFile(outfile)

    # sort lines by element and ion and energy
    lines.sort(key=lambda x: (x[0], x[1], 1./x[2]))

    # construct up FITS table to APEC format
    col_lambda = fits.Column(
        name='Lambda', format='1E', unit='A',
        array=[i[2] for i in lines])
    col_lambda_err = fits.Column(
        name='Lambda_Err', format='1E', unit='A',
        array=np.zeros( (len(lines),) ) + np.nan)
    col_epsilon = fits.Column(
        name='Epsilon', format='1E',
        unit='photons cm^3 s^-1',
        array=[i[3] for i in lines])
    col_epsilon_err = fits.Column(
        name='Epsilon_Err', format='1E',
        unit='photons cm^3 s^-1',
        array=np.zeros((len(lines),)) + np.nan)
    col_element = fits.Column(
        name='Element', format='1J',
        array=[i[0] for i in lines])
    col_ion = fits.Column(
        name='Ion', format='1J',
        array=[i[1] for i in lines])

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
    return tabhdu, len(lines), len(elements)

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
    """Interpret dumped continuum file."""

    print('Interpreting dumped continuum for T=%g' % T)

    # read in continum from each file
    continuua = {}

    allelements = ['H']+apec_elements
    for element in allelements:
        filename = 'tmp_conti_T%010f_%s.asc' % (T, element)
        energy, vals = readContinuum(filename)
        continuua[element] = vals
        deleteFile(filename)

    # subtract H continuum from each element except hydrogen
    # also divide by abundance these were generated at
    for element in apec_elements:
        continuua[element] -= continuua['H']
        continuua[element] /= continuum_mult

    # construct table
    col_element = fits.Column(
        name='Z', format='1J',
        array=[element_nums[i] for i in allelements])
    col_rmJ = fits.Column(
        name='rmJ', format='1J', array=np.zeros(len(allelements)))
    col_N_Cont = fits.Column(
        name='N_Cont', format='1J',
        array= [len(energy)]*len(allelements))
    col_E_Cont = fits.Column(
        name='E_Cont', format='%iE' % len(energy),
        unit='keV', array=np.resize(energy, (len(allelements), len(energy))))

    col_Continuum = fits.Column(
        name='Continuum', format='%iE' % len(energy),
        unit='photons cm^3 s^-1 keV^-1',
        array=[continuua[i] for i in allelements])
    col_Cont_Err = fits.Column(
        name='Cont_Err', format='%iE' % len(energy),
        unit='photons cm^3 s^-1 keV^-1',
        array=np.zeros( (len(allelements), len(energy)) ))

    # we make all the pseudo continuua zero
    col_N_Pseudo = fits.Column(
        name='N_Pseudo', format='1J',
        array=np.zeros(len(allelements)))

    # can't get below to work, unless I set array size to 2 columns
    # doesn't seem to matter, though
    col_E_Pseudo = fits.Column(
        name='E_Pseudo', format='2E',
        array=np.zeros( (len(allelements), 2) ))
    col_Pseudo = fits.Column(
        name='Pseudo', format='2E',
        array=np.zeros( (len(allelements), 2) ))
    col_Pseudo_Err = fits.Column(
        name='Pseudo_Err', format='2E',
        array=np.zeros( (len(allelements), 2) ))

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
    for i in allelements:
        totcoco += continuua[i].sum()
    h['TOT_COCO'] = totcoco

    return (tabhdu, len(allelements), len(energy)*len(allelements), 0)

def interpretAllLines():
    """Interpret spex dumped spectra."""

    # generate HDUs for each temperature
    hdus = []
    Nelement = []
    Nline = []
    for T in temperatures:
        hdu, numlines, numelements = interpretDumpedLines(T)
        Nline.append(numlines)
        Nelement.append(numelements)
        hdus.append(hdu)

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

def interpretAllContinuum():
    """Build up continuum output file."""

    # make continuum HDUs for each temperature
    hdus = []
    NElement = []
    NCont = []
    NPseudo = []
    for T in temperatures:
        hdu, numelem, numcont, numpseudo = interpretDumpedContinuum(T)
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
        col_kT, col_EDensity, col_NElement, col_NCont, col_NPseudo])
    tabhdu.name = 'PARAMETERS'

    # make output file containing the continuum
    hdulist = fits.HDUList([fits.PrimaryHDU(), tabhdu] + hdus)
    hdulist.writeto('%s_coco.fits' % outroot, overwrite=True)

def main():
    """Main routine."""
    generateOutput()
    interpretAllLines()
    interpretAllContinuum()

if __name__ == '__main__':
    main()
