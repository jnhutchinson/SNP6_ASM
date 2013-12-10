#!/usr/bin/python

# libraries
import argparse, glob, os.path, subprocess, re

# setup and parse command line arguments
parser = argparse.ArgumentParser(description="Extracting raw intensities from AffySNP6 CEL files", epilog="Make sure all your CEL file extensions are in UPPERCASE")
parser.add_argument('-c','--CELdir', type=str,  help='directory containing CEL files', required=True, action='store')
parser.add_argument('-a','--aptdir', type=str,  help='directory containing apt-cel-extract binary', required=True, action='store')
parser.add_argument('-d','--CDFfile', type=str,  help='full path to CDF file for AffySNP6 chip', required=True, action='store')

args = parser.parse_args()

# construct full path to apt-power tools binary
aptcelextract=os.path.join(args.aptdir, "apt-cel-extract" )

# switch to the directory with the CEL files
os.chdir(args.CELdir)


# process all CEL files in directory one by one

reps = {'CEL':'rawints.txt', '-':'_'}

def replace_all(text, dic):
    for i, j in dic.iteritems():
        text = text.replace(i, j)
    return text

for files in glob.glob("*.CEL"):
    	outputfile=replace_all(files, reps)
    	subprocess.call([aptcelextract, files, '-o', outputfile, '-d', args.CDFfile ])
