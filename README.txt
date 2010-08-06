- Code is organized into four directories according to simulator, 
gen = genesis
mus = moose
nrn = neuron
py = pymoose



py
 - channel.py: channel base class

 - kchans.py : all K+ channel definitions

 - nachans.py: all Na+ channel definitions

 - compartment.py: an extension of compartment class

 - cachans.py: Ca2+ channels

 - capool.py: extended version of CaConc

 - archan.py: the combined cation current

 - test.py: code for running the simulation

 - connection.py: network information in the form of graphs using
   networkX module.

nrn

 - test.hoc: code for running the simulation rest of the files are
  from original model

2009-04-27 22:03:31 (+0530)
	   Currently py/data/<date> is where the python code dumps its
	   data.
	   nrn/mydata/Vm.plot is the recorded membrane potential for
	   NEURON simulation.
