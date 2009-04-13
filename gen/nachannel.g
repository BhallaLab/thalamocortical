/**
   2008-12-22: Modified setup_napf_ss to incorporate shifts.  Added
               one shifted version of naf2 - which will be
               incorporatded into naf2 later.
 */
include utility
/********************************************************
 * fast Na channels
 ********************************************************/
function calc_naf_taum( v )
    float v
    if ( {v} < -30e-3 ) 
        return { 1.0e-3 * ( 0.025 + 0.14 * { exp { ( v  + 30e-3 ) / 10e-3 } } ) }
    else
        return { 1.0e-3 * ( 0.02 + 0.145 * { exp { ( -v - 30e-3 )/10e-3 } } ) }
    end
end

function calc_naf_tauh( v )
    float v
    return { 1e-3*(0.15+1.15/(1.0+ {exp { (v + 37e-3) / 15e-3 }})) }
end

function calc_naf_minf( v )
    float v
    return { 1.0/( 1.0 + { exp { (-v - 38e-3)/10e-3 }}) }
end

function calc_naf_hinf( v )
    float v
    return { 1.0/(1.0+ { exp { (v + 62.9e-3) / 10.7e-3 } })}
end

function setup_naf(path)
    str path
    int i
    float v = VMIN
    create_vchan {path} 3 1
    setfield {path} Ek {E_NA}
    i = 0
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} X_A->table[{i}] { calc_naf_taum { v }}
        setfield {path} X_B->table[{i}] { calc_naf_minf { v }}          
        setfield {path} Y_A->table[{i}] { calc_naf_tauh { v }}
        setfield {path} Y_B->table[{i}] { calc_naf_hinf { v }}
        v = v + dV
    end
    tweaktau {path} X
    tweaktau {path} Y
    dump_vchan_tables {path} 1 1
    // if ( {LOG_LEVEL} > 1 )
    //     echo "setup_naf: " {path} " - finished."
    // end
end
    
/********************************************************
 * fast Na channel 
 *
 * - used in Spiny Stellate Cells and ..
 ********************************************************/
function calc_naf2_taum( v )
    float v
    if ( v < -30e-3 ) 
        return { 1.0e-3 * ( 0.0125 + 0.1525 * { exp { ( v + 30e-3 ) / 10e-3 } } ) }
    else
        return { 1.0e-3 * ( 0.02 + 0.145 * { exp { ( -v - 30e-3 ) / 10e-3 } } ) }
    end
end

function calc_naf2_tauh( v )
    float v
    return { 1e-3 * (0.225 + 1.125/(1.0 + {exp { (v + 37e-3) / 15e-3 }})) }
end

// function calc_naf2_minf( v )
//     float v
//     return { 1.0 / (1.0 + {exp { (-v - 38e-3) / 10e-3 }}) }
// end

function calc_naf2_hinf(v)
    float v
    return { 1.0 / (1.0 + {exp { (v + 58.3e-3) / 6.7e-3 } })}
end

function setup_naf2(path)
    str path
    int i
    float v = VMIN
    create_vchan {path} 3 1
    setfield {path} Ek {E_NA}
    i = 0
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} X_A->table[{i}] { calc_naf2_taum { v } }
        setfield {path} X_B->table[{i}] { calc_naf_minf { v } }          
        setfield {path} Y_A->table[{i}] { calc_naf2_tauh { v } }
        setfield {path} Y_B->table[{i}] { calc_naf2_hinf { v } }
        v = v + dV
    end
    tweaktau {path} X
    tweaktau {path} Y
    dump_vchan_tables {path} 1 1
    
    if ( {LOG_LEVEL} > 1 )
        echo "setup_naf2: " {path} " - finished."
    end
end

/**
   voltage shifted NaF2
*/
function setup_naf2_shifted(path, shift)
    str path
    float shift
    int i
    float v = VMIN 
    create_vchan {path} 3 1
    setfield {path} Ek {E_NA}
    i = 0
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} X_A->table[{i}] { calc_naf2_taum { v + shift  } }
        setfield {path} X_B->table[{i}] { calc_naf_minf { v + shift } }          
        setfield {path} Y_A->table[{i}] { calc_naf2_tauh { v  } }
        setfield {path} Y_B->table[{i}] { calc_naf2_hinf { v  } }
        v = v + dV
    end
    tweaktau {path} X
    tweaktau {path} Y
    dump_vchan_tables {path} 1 1
    
    if ( {LOG_LEVEL} > 1 )
        echo "setup_naf2_shifted: " {path} " - finished."
    end
    
end

        

/********************************************************
 * fast Na channels for thalamo-cortical relay cells
* same as naf, only -3e-3 V shifted
 ********************************************************/

function setup_naf_tcr(path)
    str path
    int i
    float v = VMIN - 3e-3
    create_vchan {path} 3 1
    setfield {path} Ek {E_NA}
    i = 0
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} X_A->table[{i}] { calc_naf_taum { v }}
        setfield {path} X_B->table[{i}] { calc_naf_minf { v }}          
        setfield {path} Y_A->table[{i}] { calc_naf_tauh { v }}
        setfield {path} Y_B->table[{i}] { calc_naf_hinf { v }}
        v = v + dV
    end
    tweaktau {path} X
    tweaktau {path} Y

    // if ( {LOG_LEVEL} > 2 )
    //     dump_vchan_tables {path} 1 1
    // end

    // if ( {LOG_LEVEL} > 1 )
    //     echo "setup_naf: " {path} " - finished."
    // end
end

/****************************************************************
 * persistent Na channels for spiny thalamo-cortical relay cells
 ****************************************************************/
function calc_nap_minf(v)
    float v
    return 1.0/( 1.0 + { exp { (-v - 40e-3) / 10e-3 } } )
end

function calc_nap_taum( v )
    if ( {v} < -40e-3 ) 
        return { 1.0e-3 * ( 0.025 + 0.14 * { exp { ( v  + 40e-3 ) / 10e-3 } } ) }
    else
        return { 1.0e-3 * ( 0.02 + 0.145 * { exp { ( -v - 40e-3 )/10e-3 } } ) }
    end
end

function setup_nap(path)
    str path
    int i
    float v = VMIN
    create_vchan {path} 1
    setfield {path} Ek {E_NA}
    i = 0
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} X_A->table[{i}] { calc_nap_taum { v }}
        setfield {path} X_B->table[{i}] { calc_nap_minf { v }}          
        v = v + dV
    end
    tweaktau {path} X

    // if ( {LOG_LEVEL} > 2 )
    //     dump_vchan_tables {path} 1 0
    // end

    // if ( {LOG_LEVEL} > 1 )
    //     echo "setup_nap: " {path} " - finished."
    // end
end
/****************************************************************
 * fast persistent Na channels for spiny thalamo-cortical relay cells
 ****************************************************************/
function setup_napf(path)
    str path
    int i
    float v = VMIN
    create_vchan {path} 3
    setfield {path} Ek {E_NA}
    i = 0
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} X_A->table[{i}] { calc_naf_taum { v }}
        setfield {path} X_B->table[{i}] { calc_naf_minf { v }}          
        v = v + dV
    end
    tweaktau {path} X

    // if ( {LOG_LEVEL} > 2 )
    //     dump_vchan_tables {path} 1 0
    // end

    // if ( {LOG_LEVEL} > 1 )
    //     echo "setup_napf: " {path} " - finished."
    // end
end

/****************************************************************
 * fast persistent Na channels for spiny stellate cells
 ****************************************************************/
function setup_napf_ss(path, shift)
    str path
    float shift
    int i
    float v = VMIN + shift
    create_vchan {path} 3 0
    setfield {path} Ek {E_NA}
    i = 0
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} X_A->table[{i}] { calc_naf_taum { v } }
        setfield {path} X_B->table[{i}] { calc_naf_minf { v } }          
        v = v + dV
    end
    tweaktau {path} X

    // if ( {LOG_LEVEL} > 2 )
    //     dump_vchan_tables {path} 1 0
    // end

    // if ( {LOG_LEVEL} > 1 )
    //     echo "setup_napf: " {path} " - finished."
    // end
end

/****************************************************************
 * fast persistent Na channels for thalamo-cortical relay cells
 ****************************************************************/

function setup_napf_tcr(path)
    str path
    int i
    float v = VMIN + 7e-3
    create_vchan {path} 1
    setfield {path} Ek {E_NA}
    i = 0
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} X_A->table[{i}] { calc_naf_taum { v } }
        setfield {path} X_B->table[{i}] { calc_naf_minf { v } }          
        v = v + dV
    end
    tweaktau {path} X

    // if ( {LOG_LEVEL} > 2 )
    //     dump_vchan_tables {path} 1 0
    // end

    // if ( {LOG_LEVEL} > 1 )
    //     echo "setup_napf: " {path} " - finished."
    // end
end
