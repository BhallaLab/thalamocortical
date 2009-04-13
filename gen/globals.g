float PI = 3.141592 

str SIMULATOR

if ({version} < 3.0)
    SIMULATOR = "genesis"
else
    SIMULATOR = "moose"
end

float SIMDT = 1e-6 // !! CONSTANT !!

float PLOTDT = 1e-4 // !! CONSTANT !!

float SIMTIME // !! SET IN ONLY ONE PLACE !!

int PLOTSTEPS // !! SHOULD BE CALCULATED ONLY AFTER SIMTIME IS SET!!

int SIMSTEPS  // !! SHOULD BE CALCULATED ONLY AFTER SIMTIME IS SET!!

float VMIN = -120e-3 // xmin for creating tables

float VMAX = 40e-3  // xmax for creating tables

float NDIVS = 640 // divisions in channel gate tables

float dV = (VMAX - VMIN) / NDIVS

float INJECTION = 1e-11 // injection current - set once

float E_NA = 50e-3
float E_K = -95e-3
float E_K_FS = -100e-3 // For fast spiking cells
float E_CA = 125e-3
float E_AR = -35e-3

float TEST_RA = 2.5
float TEST_RM = 5.0
float TEST_CM = 0.009
float TEST_EM = -65e-3
float TEST_EK = -100e-3
float TEST_DIA = 15e-6
float TEST_LEN = 20e-6
// following test values are for spiny stellate level 1 compartment
// unit in Fortran code seems to be mS/cm^2
// unit in NEURON HOC code is S/cm^2
// unit here in GENESIS/MOOSE is S/m^2
float TEST_GNAF_DENS = 1500.0 
float TEST_GKDR_DENS = 1000.0
float TEST_GNAP_DENS = 0.001 * TEST_GNAF_DENS
float TEST_GCAT_DENS = 1.0
float TEST_GCAL_DENS = 5.0
float TEST_GKA_DENS = 300.0
float TEST_GKC_DENS = 100.0
float TEST_GKM_DENS = 37.5
float TEST_GK2_DENS = 1.0
float TEST_GKAHP_DENS = 1.0
float TEST_GAR_DENS = 2.5

// float PRE_TIME = 0.05
// float PULSE_WIDTH = 0.05
// float POST_TIME = 0.05
// float INJECTION = 0.3e-9

// SIMSTEPS = SIMTIME / SIMDT + 1
// PLOTSTEPS = SIMTIME / PLOTDT + 1
