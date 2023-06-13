#' default log config
#' @importFrom logger layout_glue_generator appender_file 
default_log_config <- function()
{
    layout <- layout_glue_generator("{time} | {level} | {namespace} | {msg}")
    appender <- appender_file(file = "debug.log")
    logger <- logger(appenders = appender,
        layout = layout, level = log_level_info)
    logger
}

#' @importFrom logger log_info
info <- function(logger_config = default_log_config, ...) {
  log(logger_config, log_info, ...)
}

#' @importFrom logger log_warn 
warn <- function(logger_config = default_log_config, ...) {
  log(logger_config, log_warn, ...)
}

#' @importFrom logger log_error
error <- function(logger_config = default_log_config, ...) {
  log(logger_config, log_error, ...)
}

log <- function(log_func, logger_config = default_log_config(), ...) {
  logger <- logger_config()
  
  send_message <- function(x, msg, log_func) {
    if (!is.null(x)) {
      if (is.integer(x) || is.numeric(x)) {
        log_func(logger, sprintf("| %s: %s", msg, toString(x)))
      } else if (is.matrix(x) || is.vector(x) || is.data.frame(x)) {
        log_func(logger, sprintf("| %s: %s", msg, toString(head(x))))
      } else if (is.list(x)) {
        log_func(logger, sprintf("| %s is a list of length %d, with element classes: %s", msg, length(x), toString(sapply(x, class))))
      } else {
        log_func(logger, sprintf("| %s: %s", msg, toString(x)))
      }
    }
  }
  
  items <- list(...)
  item_names <- names(items)
  for (i in 1:length(items))
    send_message(items[[i]], item_names[i], log_func)
}
