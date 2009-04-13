/*******************************************************************
 * File:            globals.g
 * Description:     Contains the global variable list
 *                  This file must be included before anything else.
 * Author:          Subhasis Ray
 * E-mail:          ray dot subhasis at gmail dot com
 * Created:         2008-11-29 02:14:18 (+0530)
 ********************************************************************/

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
