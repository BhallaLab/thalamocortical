/*******************************************************************
 * File:            kchannel.g
 * Description:      
 * Author:          Subhasis Ray
 * E-mail:          ray dot subhasis at gmail dot com
 * Created:         2008-12-05 15:36:50 (+0530)
 * 
 * Change log:
 *
 * 2008-12-05 15:37:30 (+0530)
 *      Added new vesion of create_vchan in the setup methods. Thus 
 *      the channel gets created at setup time. also the tables in
 *      gates also get allocated.
 *      Added kc_fast dynamics.
 ********************************************************************/

include utility
/**********************************************************************
 * K+ channel - delayed rectifier.  The definition for the generic
 * version is different in neuron model ( kdr.mod ) and traub's
 * fortran model (most of the integrate_XXX.f where XXX is the cell
 * type )
 **********************************************************************/

function calc_taum_kdr( v )
    float v
    if ( {v} <= -10.0e-3 )
        return 1.0e-3 * (0.25 + 4.35 * { exp { (v + 10.0e-3) / 10.0e-3 } })
    else 
        return 1.0e-3 * (0.25 + 4.35 * { exp { (-v - 10.0e-3) / 10.0e-3 } } )
    end
end

function calc_minf_kdr( v )
    float v
    return  1.0/(1.0+ { exp { (-v-29.5e-3)/10.0e-3} } )
end

function setup_kdr( path )
    str path
    float v = VMIN
    int i
    create_vchan {path} 4 0
    setfield {path} Ek {E_K}
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} \
            X_A->table[{i}] { calc_taum_kdr { v }  } \
            X_B->table[{i}] { calc_minf_kdr { v }  }          
        v = v + dV
    end
    tweaktau {path} X
    if ( {LOG_LEVEL} <= {DEBUG} )
        dump_vchan_tables {path}  1 0
    end
    // if ( {LOG_LEVEL} > 1 )
    //     echo "setup_kdr: " {path} " - finished."
    // end
end    

function calc_minf_kdr_fs( v )
    float v
    return  1.0 / (1.0 + {exp {(-v - 27e-3) / 11.5e-3 }})     
end

function setup_kdr_fs( path )
    str path
    float v = VMIN
    int i
    create_vchan {path} 4 0
    setfield {path} Ek {E_K_FS}
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} \
            X_A->table[{i}] { calc_taum_kdr { v } } \
            X_B->table[{i}] { calc_minf_kdr_fs { v } }          
        v = v + dV
    end
    tweaktau {path} X
    if ( {LOG_LEVEL} <= {DEBUG} )
        dump_vchan_tables {path}  1 0
    end
end    

/**
   This is identical to kdr - hence commenting out

function calc_minf_kdrtcr( v )
    float v
    return  1.0/(1.0+ { exp { (-v-29.5e-3)/10.0e-3 } } )     
end

function setup_kdrtcr( path )
    str path
    float v = VMIN
    int i
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} \
            X_A->table[{i}] { calc_taum_kdr { v } } \
            X_B->table[{i}] { calc_minf_kdrtcr { v } }          
        v = v + dV
    end
    tweaktau {path} X
    if ( {LOG_LEVEL} > 2 )
        dump_vchan_tables {path}  1 0
    end
    if ( {LOG_LEVEL} > 1 )
        echo "setup_kdrtcr: " {path} " - finished."
    end
end    
*/

/***************************************************** 
 * KA channel  - identical for all cells
 *****************************************************/
function calc_minf_ka( v )
    float v
    return 1.0 / (1.0+ {exp { (-v - 60.0e-3) / 8.5e-3 } })
end

function calc_taum_ka( v )
    float v
    return 0.185e-3 + 0.5e-3 / ({exp {(v + 35.8e-3) / 19.7e-3}} + {exp {(-v - 79.7e-3) / 12.7e-3}})
end

function calc_hinf_ka( v )
    float v
    return 1.0/( 1.0 + { exp {(v + 78.0e-3) / 6.0e-3}})
end

function calc_tauh_ka( v )
    float v
    if ( {v} < -63.0e-3 ) 
        return 0.5e-3/ ( { exp { (v + 46.0e-3) / 5.0e-3 } } + { exp {(-v - 238.0e-3) / 37.5e-3 } } )
    else
        return 9.5e-3
    end
end

function setup_ka( path )
    str path
    float v = VMIN
    int i
    create_vchan {path} 4 1
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} \
            X_A->table[{i}] { calc_taum_ka {v} } \
            X_B->table[{i}] { calc_minf_ka {v} } \
            Y_A->table[{i}] { calc_tauh_ka {v} } \
            Y_B->table[{i}] { calc_hinf_ka {v} }
        v = v + dV
    end
    
    tweaktau {path} X
    tweaktau {path} Y
    if ( {LOG_LEVEL} <= {DEBUG} )
        dump_vchan_tables {path}  1 1
    end
end


/**
   for tufted intrinsic bursting cells the tau_h parameter gets multiplied by 2.6
*/
function setup_ka_ib( path )
    str path
    float v = VMIN
    int i
    create_vchan {path} 4 1
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} \
            X_A->table[{i}] { calc_taum_ka {v} } \
            X_B->table[{i}] { calc_minf_ka {v} } \
            Y_A->table[{i}] { {calc_tauh_ka {v}} * 2.6 } \
            Y_B->table[{i}] { calc_hinf_ka {v} }
        v = v + dV
    end
    
    tweaktau {path} X
    tweaktau {path} Y
    if ( {LOG_LEVEL} <= {DEBUG} )
        dump_vchan_tables {path}  1 1
    end
end

/***************************************************** 
 * K2 channel 
 *****************************************************/
function calc_minf_k2 ( v )
    float v
    return 1.0/(1.0+ { exp { (-v  - 10.0e-3 )/17.0e-3 } } )
end

function calc_taum_k2( v )
    float v
    return 4.95e-3 + 0.5e-3/( {exp { ( v  - 81e-3) / 25.6e-3 }} + { exp { ( - v - 132e-3) / 18e-3 } } )
end

function calc_hinf_k2( v )
    float v
    return 1.0/( 1.0 + { exp { ( v  + 58e-3 )/10.6e-3 } } )
end

function calc_tauh_k2(v)
    float  v
    return 60e-3 + 0.5e-3/( { exp { (v - 1.33e-3)/200e-3 }} + { exp { (-v -130e-3)/7.1e-3 } } )
end


function setup_k2 ( path )
    str path
    float v = VMIN
    int i
    // if ( {LOG_LEVEL} > 1 )
    //     echo "setup_k2: setting up " {path}
    // end
    create_vchan {path} 1 1
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} \
            X_A->table[{i}] { calc_taum_k2 {v} } \
            X_B->table[{i}] { calc_minf_k2 {v} } \
            Y_A->table[{i}] { calc_tauh_k2 {v} } \
            Y_B->table[{i}] { calc_hinf_k2 {v} }
        v = v + dV
    end
    
    tweaktau {path} X
    tweaktau {path} Y
    if ( {LOG_LEVEL} <= {DEBUG} )
        dump_vchan_tables {path}  1 1
    end
    // if ( {LOG_LEVEL} > 1 )
    //     echo "setup_k2: " {path} " - finished."
    // end
end

/*********************************************************************
 * KC channel - dependent on voltage & [Ca+2 ] 
 *********************************************************************/
function setup_kc ( path )
    str path
    float v, alpham, betam, chi, dchi
    int i
    // if ( {LOG_LEVEL} > 1 )
    //     echo "setup_kc: setting up " {path}
    // end
    create_vchan {path} 1 0
    setfield {path} Zpower 1
    v = VMIN
    for ( i = 0; i <= NDIVS; i = i + 1 )
        if ( {v } < -10e-3 )
            alpham = (2e3 / 37.95) * {exp {(v + 50e-3) / 11e-3 - (v + 53.5e-3) / 27e-3}}
            betam = 2e3 * {exp {(-v - 53.5e-3) / 27e-3}} - alpham
        else
            alpham = 2e3 * {exp {(-v - 53.5e-3) / 27e-3}}
            betam = 0.0
        end
        setfield {path} \
            X_A->table[{i}] { alpham } \
            X_B->table[{i}] { betam } 
        v = v + dV
    end
    tweakalpha {path} X
    setfield {path} X_A->calc_mode 0 X_B->calc_mode 0
    // setup z-gate -- ca2+ dependency
    float xmin = 0.0
    float xmax = 250.0
    float y
    int xdivs = 3000

    setfield {path} Zpower 1 //instant {INSTANTZ}
    call {path} TABCREATE Z {xdivs} {xmin} {xmax}
    chi = xmin
    dchi = (xmax - xmin) / xdivs
    i = 0
    // bad loop ( if inside for )
    // and i am not sure about the units
    for ( i = 0; i <= {xdivs}; i = i + 1 )
        if (chi < 250.0)
            y = chi/250.0
        else
            y = 1.0
        end
        setfield {path} Z_A->table[{i}] {y}
        setfield {path} Z_B->table[{i}] 1.0
        chi = chi + dchi
    end
    setfield {path} Z_A->calc_mode 0 Z_B->calc_mode 0
    setfield {path} instant {INSTANTZ}
    call {path} TABFILL Z 3000 2
    if ( {LOG_LEVEL} <= {DEBUG} )
        dump_vchan_tables {path}  1 0
        tab2file {SIMULATOR}_kc_za.plot {path} Z_A -overwrite
        tab2file {SIMULATOR}_kc_zb.plot {path} Z_B -overwrite
        echo setup_kc: {path} - dumped Z tables in {SIMULATOR}_za_kc.plot and {SIMULATOR}_zb_kc.plot 
    end
end


/**
 kc channel with C kinetics speedup ( alpha and beta twice that of ordinary kc channel )
*/
function setup_kc_fast(path)
    str path
    int i
    setup_kc {path}
    for ( i = 0; i <= 3000; i = i + 1 )
        setfield {path} \
            X_A->table[{i}] 2 * { getfield {path} X_A->table[{i}] } \
            X_B->table[{i}] 2* { getfield {path} X_A->table[{i}] } 
    end
end


/*********************************************************************
 * KAHP channel - only [Ca+2 ] dependent 
 *********************************************************************/
function setup_kahp ( path )
    str pathtr
    float chi, dchi
    int i
    create_vchan {path} 0 0
    setfield {path} Zpower 1
    float xmin = 0.0
    float xmax = 1000.0
    int xdivs = 3000
    call {path} TABCREATE Z {xdivs} {xmin} {xmax}
    chi = 0.0
    dchi = (xmax - xmin)/xdivs
    i = 0
    // taken from neuron model kahp.mod - does not look like same as traub's:
    // alpham_ahp = min( 0.2e-4*chi, 0.01)
    // which translates to:
    // if chi < 0.5e-3 then alpha = 0.2e-4*chi; else alpha = 0.01
    // alpha is multiplied by 10^3 to convert units (ms -> s)
    for ( i = 0; i <= {xdivs}; i = i + 1 )
        if ( {chi} < 100.0 )
            setfield {path} \
                Z_A->table[{i}] {0.1 * chi} 
        else
            setfield {path} \
                Z_A->table[{i}] 10.0 
        end // if
        setfield {path} Z_B->table[{i}] 10.0
        chi = chi + dchi
    end // for
    setfield {path} Z_A->calc_mode 0 Z_B->calc_mode 0
    tweakalpha {path} Z
    if ( {LOG_LEVEL} <= {DEBUG} )
        tab2file {SIMULATOR}_za_kahp.plot {path} Z_A -overwrite 
        tab2file {SIMULATOR}_zb_kahp.plot {path} Z_B -overwrite
        echo setup_kahp: {path} - dumped Z tables in {SIMULATOR}_za_kahp.plot and {SIMULATOR}_zb_kc.plot
    end
end


/*********************************************************************
 * KAHP slower channel - only [Ca+2 ] dependent 
 *********************************************************************/
function setup_kahp_slow ( path )
    str path
    float chi, dchi
    int i
    create_vchan {path} 0 0
    setfield {path} Zpower 1
    float xmin = 0.0
    float xmax = 500.0
    int xdivs = 3000
    call {path} TABCREATE Z {xdivs} {xmin} {xmax}
    chi = 0.0
    dchi = (xmax - xmin)/xdivs
    i = 0
    // taken from neuron model kahp_slower.mod
    // alpham_ahp = min( 0.2e-4*chi, 0.01)
    // which translates to:
    // if chi < 0.5e-3 then alpha = 0.2e-4*chi; else alpha = 0.01
    // alpha is multiplied by 10^3 to convert units (ms -> s)
    for ( i = 0; i <= {xdivs}; i = i + 1 )
        if ( {chi} < 500.0 )
            setfield {path} \
                Z_A->table[{i}] {0.02*chi} 
        else
            setfield {path} \
                Z_A->table[{i}] 10.0 
        end // if
        setfield {path} Z_B->table[{i}] 10.0
        chi = chi + dchi
    end // for
    setfield {path} Z_A->calc_mode 0 Z_B->calc_mode 0
    tweakalpha {path} Z
    if ( {LOG_LEVEL} <= {DEBUG} )
        tab2file {SIMULATOR}_za_kahp_slow.plot {path} Z_A -overwrite 
        tab2file {SIMULATOR}_zb_kahp_slow.plot {path} Z_B -overwrite
        echo setup_kahp_slow: {path} - dumped Z tables in {SIMULATOR}_za_kahp.plot and {SIMULATOR}_zb_kc.plot
    end
end


/*********************************************************************
 * KM channel - common for TCR, SS
 *********************************************************************/
function calc_alpham_km ( v )
    float v
    return 0.02 / ( 1 + { exp {  ( -v - 20e-3 ) / 5e-3 } } )
end

function calc_betam_km( v )
    float v
    return 0.01 * { exp { ( -v - 43e-3 ) / 18e-3 } }
end


function setup_km ( path )
    str path
    float v = VMIN
    int i
    create_vchan {path} 1 0
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} \
            X_A->table[{i}] { calc_alpham_km {v} } \
            X_B->table[{i}] { calc_betam_km {v} } 
        v = v + dV
    end
    tweakalpha {path} X
    if ( {LOG_LEVEL} <= {DEBUG} )
        dump_vchan_tables {path}  1 0
    end
end


/**
   ar channel - the combined cation channel 
*/
function setup_ar(path)
    str path
    float v = VMIN
    int i
    create_vchan {path} 1 0
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} \
            X_A->table[{i}] { 1.0e-3 / { exp {( v + 75e-3) / 5.5e-3 } } } \
            X_B->table[{i}] { 1.0e-3 / {{ exp { -14.6 - 86 * v } } + { exp { - 1.87 + 70 * v } }}}
        v = v + dV
    end
    tweaktau {path} X
    if ( {LOG_LEVEL} <= {DEBUG} )
        dump_vchan_tables {path} 1 0
    end
end
