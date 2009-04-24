- Code is organized into four directories according to simulator, 
gen = genesis
mus = moose
nrn = neuron
py = pymoose

each directory has a 'data' subdirectory.

py
 - kchans.py : all K+ channel definitions
 - nachans.py: all Na+ channel definitions
 - compartment.py: an extension of compartment class
 - cachans.py: Ca2+ channels
 - capool.py: extended version of CaConc
 - test.py: code for running the simulation

nrn
 - test.hoc: code for running the sunulation 
 rest of the files are from original model
