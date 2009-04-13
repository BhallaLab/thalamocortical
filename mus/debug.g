/*******************************************************************
 * File:            debug.g
 * Description:     functions for debugging. mainly log various 
 *                  parameters.
 * Author:          Subhasis Ray
 * E-mail:          ray dot subhasis at gmail dot com
 * Created:         2008-11-29 02:30:14 (+0530)
 ********************************************************************/

include globals
include logging

// list of relevant fields of a compartment
str COMPARTMENT_FIELDS = "len dia initVm Vm Em Rm Cm Ra inject Im" 

// list of relevant fields of a channel
str CHANNEL_FIELDS = "Ek Gbar Xpower Ypower Zpower instant Ik"

// list of relevant fields of a pulsegen
str PULSEGEN_FIELDS = "level1 width1 delay1 level2 width2 delay2 baselevel trig_time trig_mode"

/**
   log the fields of an object.
*/
function log_fields(object, field_list)
    str object // path of the object to be shown
    str field_list // a string containing the list of fields to be logged

    str field
    foreach field ( {arglist {field_list} } )
        log_info {{object} @ "." @ {field} @ " = " @ {getfield {object} {field}}}
    end
end

/**
   save data tables within a container to files.
*/
function dump_tables(data_container)
    str data_container // the element containing the table elements

    str file_name
    str data_table
    foreach data_table ({el {data_container}/#})
        file_name = {SIMULATOR} @ "_" @ {getfield {data_table} name} @ ".plot"
        tab2file {file_name} {data_table} -overwrite
    end
end
