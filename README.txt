This is the PyMOOSE version of Single-column thalamocortical model by
Traub et al 2005 with enhancements to add variability in single cell
properties as well as synaptic condctances. 

Code is organized into four directories according to language/simulation system:

for: original fortran code

g77traub: My update of the g77 port to compile and run on mpif77.

nrn: neuron model with additional scripts/updates for dumping information and testing.

py: pymoose version. This evolved with moose from 2009 till 2013 when a
    drastic overhaul in te API was made (dh_branch). By the time the
    new stuff stabilized, it was too much work to convert the entire
    model into this new API. Still, the single cell models and ion
    channels have been translated and are available as part of
    moose-examples. They are yet to be extensively tested.

    The files in this diectory are as follows:

    {Cell}.p - prototype file for cells of type {Cell}. These are in
    	       GENESIS prototype format and are read by moose to
    	       create the prototype cell model. The prototype is
    	       copied in moose to make the cell population.  py -
    	       channel.py: channel base class

    {cell}.py - Script defining the Python class for Cell. The
	        __init__ method can take an existing object and a path
	        and make a copy of the existing object in the
	        specified path.

    trbsim.py - main script for running the simulation. It takes a
		large number of options.

    trbnetdata.py: network information from Traub et al 2005.

    trbnet.py - class TraubNet in this script creates the actual
    	      	thalamocortical network based on specified
    	      	cell-to-cell or celltype-to-celltype graphs or with
    	      	default values from the original Traub et al 2005
    	      	model.

    connection.py - superseeded by trbnet.py			

    kchans.py : all K+ channel definitions
    
    nachans.py: all Na+ channel definitions
    
    compartment.py: an extension of compartment class
    
    cachans.py: Ca2+ channels
    
    capool.py: extended version of CaConc
    
    archan.py: the combined cation current

    defaults.ini: default simulation configuration parameters.

    custom.ini: customizations to override the default simulation
    		parameters.

Currently py/data/<date> is where the python code dumps its data.
