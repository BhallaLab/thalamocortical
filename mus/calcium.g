/*******************************************************************
 * File:            calcium.g
 * Description:     Ca2+ ion channels and Ca2+ concentration
 * Author:          Subhasis Ray
 * E-mail:          ray dot subhasis at gmail dot com
 * Created:         2008-12-11 20:25:42 (+0530)
 ********************************************************************/

/**
   create and setup a Ca_conc object
*/
function setup_ca_conc(path, thickness, tau)
    str path // path of the Ca_conc object to be created
    float thickness // thickness of the shell near the membrane surface
    float tau // time constant for Ca2+ decay
    float dia, len

    create Ca_concen {path}
    pushe {path}
    dia = {getfield .. dia}
    len = {getfield .. len}
    pope
    setfield {path} \
        tau {tau} \
        B {1 / (2 * 96485 * PI * dia * len * thickness)} \
        Ca_base 0.0
end


/**
   cal - long lasting Ca2+ currents.
*/
function calc_alpham_cal( v )
    float v
    return 1.6e3/( 1.0+ { exp { -0.072e3 * ( v - 5.0e-3 ) } } )
end

function calc_betam_cal( v )
    float v, tmp
    tmp = { v } + 8.9e-3
    if ({abs {tmp}} < 1e-9)
        return  0.1e3 * {exp {- {tmp} / 5.0}}
    else
        return 0.02e3 * {tmp} / ({exp {{tmp} / 5.0}} - 1.0)
    end
end

function setup_cal( path )
    str path
    float v = VMIN
    int i

    create_vchan {path} 2 0
    for ( i = 0; i <= NDIVS; i = i + 1 )
        setfield {path} \
            X_A->table[{i}] { calc_alpham_cal { v }  } \
            X_B->table[{i}] { calc_betam_cal { v }  }          
        v = v + dV
    end
    tweakalpha {path} X
    if ({LOG_LEVEL} <= {DEBUG})
        dump_vchan_tables {path} 1 0
    end
end    


/***************************************************** 
 * CAT channel 
 *****************************************************/
function calc_minf_cat ( v )
    float v
    return 1.0/(1.0+ { exp { (-v-56.0e-3)/6.2e-3 } } )
end

function calc_taum_cat( v )
    float v
    return 1.0e-3*(0.204 + 0.333/( { exp {(v+15.8e-3)/18.2e-3} } + { exp { (-v-131.0e-3)/16.7e-3 }}))
end

function calc_hinf_cat( v )
    float v
    return 1.0/(1.0+ { exp { (v+80.0e-3)/4.0e-3 } } )
end

function calc_tauh_cat(v)
    float  v
    if ( {v} < -81.0e-3 )
       return 1e-3 * 0.333 * { exp {(v + 466e-3) / 66.6e-3} }
    else
       return 1e-3 * (9.32 + 0.333 * { exp { (-v-21.0e-3) / 10.5e-3 } } )
    end
end

function setup_cat ( path )
    str path
    float v = VMIN
    int i
    create_vchan {path} 2 1
    for (i = 0; i <= NDIVS; i = i + 1)
        setfield {path} \
            X_A->table[{i}] { calc_taum_cat {v} } \
            X_B->table[{i}] { calc_minf_cat {v} } \
            Y_A->table[{i}] { calc_tauh_cat {v} } \
            Y_B->table[{i}] { calc_hinf_cat {v} }
        v = v + dV
    end
    
    tweaktau {path} X
    tweaktau {path} Y
    if ( {LOG_LEVEL} <= {DEBUG} )
        dump_vchan_tables {path} 1 1
    end

end


/**
   CaT channel variation
*/
function setup_cat_a(path)
    str path
    create_vchan {path} 2 1
    
    float v = VMIN
    int i
    for (i = 0; i <= NDIVS; i = i + 1)
        setfield {path} \
            X_A->table[{i}] {1e-3 * (1 + 0.33 / ({exp {({v} + 27e-3) / 10e-3}} + {exp {(-{v} - 102e-3) / 15e-3}}))} \
            X_B->table[{i}] {1.0 / (1.0 + {exp {(-{v} - 52e-3) / 7.4e-3}})} \
            Y_A->table[{i}] {1e-3 * (28.3 + 0.33 / ({exp {({v} + 48e-3) / 4e-3}} + {exp {(-{v} - 407e-3) / 50e-3}}))} \
            Y_B->table[{i}] {1.0 / (1 + {exp {({v} + 80e-3) / 5e-3}})}
         v = v + dV
    end
    tweaktau {path} X
    tweaktau {path} Y
    if ({LOG_LEVEL} <= {DEBUG})
        dump_vchan_tables {path} 1 1
    end
end
