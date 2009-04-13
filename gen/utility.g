/*******************************************************************
 * File:            utility.g
 * Description:     
 * Author:          Subhasis Ray
 * E-mail:          ray dot subhasis at gmail dot com
 * Created:         2008-11-29 04:32:08 (+0530)
 * 
 * ChangeLog:
 * 2008-12-05 14:53:31 (+0530) 
 * create_vchan: Removed the parameter Ek. Now it creates channel, 
 *              sets Xpower, Ypower and calls TABCREATE. This 
 *              function will be called inside setup_XYZ where XYZ
 *              is a channel. Thus the Xpower and Ypower information
 *              remains encapsulated in the channel creation and 
*               setup function.
 * 
 ********************************************************************/

include logging

/**
   Insert a current injection into a specified compartment. The object
   is created under the compartment and is named injection.

   comp_path - path of the compartment into which injection is to be
               applied.
               
   delay (s) - delay from start of simulation till the start of injection

   width (s) - duration of injection 

   level (A) - amplitude of injection current 
*/
function insert_injection(comp_path, delay, width, level)
    str comp_path
    float delay, width, level

    create pulsegen {comp_path}/injection
    setfield ^ \
        delay1 {delay} \
        width1 {width} \
        level1 {level}

    addmsg ^ {comp_path} INJECT output
end


/**
   Create a channel with reversal potential Ek
*/
function create_vchan(path, xpower, ypower)
    str compartment, path
    int xpower, ypower
    create tabchannel {path}
    setfield ^  Xpower {xpower} Ypower {ypower}
    if ( {xpower} != 0 )
        call {path} TABCREATE X {NDIVS} {VMIN} {VMAX}
    end
    if ( {ypower} != 0 )
        call {path} TABCREATE Y {NDIVS} {VMIN} {VMAX}
    end
end        


/**
   set the Gbar for channel from conductance density gbar_dens
*/
function set_gbar(channel, gbar_dens)
    str channel
    float gbar_dens
    float len, dia

    pushe {channel} // need this to bypass MOOSE incompatibility - .. not allowed except for the first item in path
    len = {getfield .. len}
    dia = {getfield .. dia}
    setfield . Gbar {gbar_dens * PI * dia * len}
    str msg = {getfield . name} @ ".Gbar set to" @ {getfield . Gbar}
    echo {msg}
    pope
end


/**
   sets up hhchannels inside the compartment.
   assumes Gbar of the channels is the specific value and turns it into absolute value.
   connects the channel to the compartment
*/
function setup_hhchannels(comp)
    str comp
    str chan
    foreach chan ({el {comp}/##[TYPE=tabchannel]})
        setfield {chan} Gbar {{getfield {chan} Gbar} * PI * {getfield {comp} dia} * {getfield {comp} len}}
        addmsg {chan} {comp} CHANNEL Gk Ek
        addmsg {comp} {chan} VOLTAGE Vm
    end
end
  
  
/**
   dump the xGate and yGate tables in channel
*/
function dump_vchan_tables(channel, x, y)
    str channel
    int x, y
    str prefix = {SIMULATOR} @ "_" @ {getfield {channel} name}
    echo dump_vchan_tables: start
    
    if ( {x} != 0 )
        tab2file {prefix}_xa.plot {channel} X_A -overwrite 
        tab2file {prefix}_xb.plot {channel} X_B -overwrite
        // if ( {LOG_LEVEL} > 1 )
        //     echo dump_vchan_tables: {channel} - dumped X tables in {SIMULATOR}_xa_{getfield {channel} name}.plot and {SIMULATOR}_xb_{getfield {channel} name}.plot
        // end
    end
    if ( {y} != 0 )
        tab2file {prefix}_ya.plot {channel} Y_A -overwrite
        tab2file {prefix}_yb.plot {channel} Y_B -overwrite        
        // if ( {LOG_LEVEL} > 1 )
        //     echo dump_vchan_tables: {channel} - dumped Y tables in {SIMULATOR}_ya_{getfield {channel} name}.plot and {SIMULATOR}_yb_{getfield {channel} name}.plot
        // end
    end
    echo dump_vchan_tables: end
end


/**
   create tables for recording channel currents
*/
function create_channel_output(channel, data)
        str channel, data
        str prefix
        pushe {channel}
        prefix = {getfield .. name} @ "_" @ { getfield . name }
        pope
        create table {data}/{prefix}_Ik
        setfield ^ step_mode 3
        call ^ TABCREATE {PLOTSTEPS} 0.0 1.0
        addmsg {channel} ^ INPUT Ik
        // if ( {VERBOSITY} > 1 )
        //     echo "create_chan_output: Created table: " /data/I_{suffix} to read {channel}.Ik
        // end
 
        create table {data}/{prefix}_Gk
        setfield ^ step_mode 3
        call ^ TABCREATE {PLOTSTEPS} 0.0 1.0
        addmsg {channel} ^ INPUT Gk
        // if ( {VERBOSITY} > 1 )
        //     echo "create_chan_output: Created table: " /data/G_{suffix} to read {channel}.Gk
        // end
        // if ( {VERBOSITY} > 3 )
        //     echo create_chan_output: END
        // end
end

function create_electronics(path)
        str path
        if(!{exists {path}})
            create neutral {path}
        end
        echo "Creating electronics."
        pushe {path}
        create pulsegen pulsegen
        echo "Created pulsegen:" {path}/pulsegen
        create diffamp iclamp
        setfield ^ gain 0.0 saturation 999.0
        echo "Created iclamp:" {path}/iclamp "with gain="{getfield iclamp gain}  "saturation="{getfield iclamp saturation}
        create diffamp vclamp
        setfield ^ gain 0.0 saturation 999.0
        echo "Created vclamp:" {path}/vclamp "with gain="{getfield vclamp gain}  "saturation="{getfield vclamp saturation}
        create RC lowpass
        setfield ^ R 1.0 C 1e-6
        echo "Created lowpass filter:" {path}/lowpass "with R="{getfield lowpass R}  "C="{getfield lowpass C}
        create PID PID
        setfield ^ gain 0.5 tau_i 20e-6 tau_d 5e-6
        echo "Created PID:" {path}/PID "with gain="{getfield PID gain} "tau_i="{getfield PID tau_i} "tau_d="{getfield PID tau_d}
        addmsg pulsegen iclamp PLUS output
        echo "Connected pulsegen to iclamp."
        addmsg pulsegen lowpass INJECT output
        echo "Connected pulsegen to lowpass"
        addmsg lowpass vclamp PLUS state
        echo "Connected lowpass to vclamp"
        addmsg vclamp PID CMD output
        echo "Connected vclamp to PID"
        pope
end

function connect_electronics(elec_path, comp_path)
        str elec_path // path of the electronics container
        str comp_path // path of the compartment to be clamped
        addmsg {elec_path}/PID {comp_path} INJECT output
        echo "Connected PID/output to "{comp_path}"/injectDest"
        addmsg {comp_path} {elec_path}/PID SNS Vm
        echo "Connected" {comp_path}"/Vm to PID/SNS"
        addmsg {elec_path}/iclamp {comp_path} INJECT output
        echo "Connected" {elec_path}"/iclamp to" {comp_path}"/injectDest"
end

function set_vclamp(elec_path, pre_t, pre_v, clamp_t, clamp_v)
        str elec_path
        float pre_v, pre_t, clamp_v, clamp_t
        setfield {elec_path}/iclamp gain 0.0
        setfield {elec_path}/vclamp gain 1.0 
        setfield {elec_path}/PID gain 0.5
        setfield {elec_path}/pulsegen baselevel {pre_v} delay1 {pre_t} level1 {clamp_v} width1 {clamp_t} delay2 1e6
        echo "Set vclamp with pre_pulse voltage:"{getfield {elec_path}/pulsegen baselevel} "pre_pulse time:"{getfield {elec_path}/pulsegen delay1} "clamp voltage:"{getfield {elec_path}/pulsegen level1} "clamp time: "{getfield {elec_path}/pulsegen width1}
end

function set_iclamp(elec_path, pre_t, inject, dur)
        str elec_path
        float pre_t, inject, dur
        setfield {elec_path}/iclamp gain 1.0
        setfield {elec_path}/vclamp gain 0.0
        setfield {elec_path}/PID gain 0.0
        setfield {elec_path}/pulsegen baselevel 0 delay1 {pre_t} level1 {inject} width1 {dur} delay2 1e6
end
