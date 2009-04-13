/*******************************************************************
 * File:            logging.g
 * Description:     function for logging to a file
 * Author:          Subhasis Ray
 * E-mail:          ray dot subhasis at gmail dot com
 * Created:         2008-11-17 20:13:03 (+0530)
 ********************************************************************/


int CRITICAL = 50
int ERROR =	40
int WARNING = 30
int INFO = 20
int DEBUG = 10
int NOTSET = 0
int LOG_LEVEL = WARNING
str LOG_FILE_PATH = "stdout"

/**
   setup the log file and log level.  must be called before using any
   logging functions.
 */
function log_init(level, file_path, mode)
    str file_name
    int level
    str mode

    LOG_LEVEL = level
    if ({LOG_FILE_PATH} != "stdout" && LOG_FILE_PATH != "")
        closefile {LOG_FILE_PATH}
    end
    LOG_FILE_PATH = {file_path}
    openfile {LOG_FILE_PATH} {mode}
end

/**
   write a log entry.  is used in other log functions for level
   specific logging.
*/
function log(level, message)
    int level
    str message
        
    if ({level} >= {LOG_LEVEL})
        if ({LOG_FILE_PATH} != "stdout" )
            writefile {LOG_FILE_PATH} {message}
        else
            echo {message}
        end
    end
end

/**
   write a debug log
*/
function log_debug(message)
    str message
    str new_message = "DEBUG " @ {message} 
    log {DEBUG}  {new_message}
end
        
/**
   write an info log
*/
function log_info(message)
    str message
    str new_message = "INFO " @ {message}
    log {INFO} {new_message}
end

/**
   write a warn log
*/
function log_warning(message)
    str message
    str new_message = "WARNING " @ {message}
    log {WARNING} {new_message}
end

/**
   write a error log
*/
function log_error(message)
    str message
    str new_message = "ERROR " @ {message}
    log {ERROR} {new_message}
end

/**
   write a critical log
*/
function log_critical(message)
    str message
    str new_message = "CRITICAL " @ {message}
    log {CRITICAL} {new_message}
end

/**
   close the log file should be called before exit. Though by default
   MOOSE and GENESIS close all open files.
 */
function log_end
    if ({LOG_FILE_PATH} != "stdout")
        closefile {LOG_FILE_PATH}
        LOG_FILE_PATH = "stdout"
    end
end


function test_log
    log_init {DEBUG} "test.log" "w"
    log_debug "Hello"
    log_info "World"
    log_warning "War"
    log_error "of the worlds"
    log_critical "is over."
    log_end
end
    
