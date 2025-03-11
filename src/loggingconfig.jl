module LoggingConfig

using Logging, LoggingExtras, Crayons

# Custom log formatting function
function custom_log_formatter(io, args)
    level = args.level
    message = args.message
    file = args.file
    line = args.line

    color = if level == Logging.Info
        Crayon(foreground=:green)
    elseif level == Logging.Warn
        Crayon(foreground=:yellow)
    elseif level == Logging.Error
        Crayon(foreground=:red)
    elseif level == Logging.Debug
        Crayon(foreground=:cyan)
    else
        Crayon(foreground=:white)
    end

    println(io, color, "[", level, "] ", message, " (", file, ":", line, ")", Crayon(reset=true))
end

# Setup custom logger
function setup_logger()
    global_logger(FormatLogger(custom_log_formatter))
end

export setup_logger

end # module LoggingConfig
