#include "Logger.h"

String LEVEL_TRACE = "TRACE";
String LEVEL_DEBUG = "DEBUG";
String LEVEL_INFO = "INFO";
String LEVEL_SUCCESS = "SUCCESS";
String LEVEL_WARNING = "WARNING";
String LEVEL_ERROR = "ERROR";
String LEVEL_CRITICAL = "CRITICAL";

String LOG_LEVEL_STRINGS[] = {LEVEL_TRACE,
                              LEVEL_DEBUG,
                              LEVEL_INFO,
                              LEVEL_SUCCESS,
                              LEVEL_WARNING,
                              LEVEL_ERROR,
                              LEVEL_CRITICAL};

Logger::Logger() : _level(SUCCESS)
{
}

void Logger::begin(uint32_t baud)
{
    _serial->begin(baud);
}

void Logger::setSerial(HardwareSerial &serial)
{
    _serial = &serial;
}

HardwareSerial &Logger::getSerial()
{
    return *_serial;
}

void Logger::setLogLevel(Level level)
{
    _level = level;
}

Level Logger::getLogLevel()
{
    return _level;
}

void Logger::trace(const char *message)
{
    log(TRACE, message);
}

void Logger::trace(const String message)
{
    log(TRACE, message);
}

void Logger::debug(const char *message)
{
    log(DEBUG, message);
}

void Logger::debug(const String message)
{
    log(DEBUG, message);
}

void Logger::info(const char *message)
{
    log(INFO, message);
}

void Logger::info(const String message)
{
    log(INFO, message);
}

void Logger::success(const char *message)
{
    log(SUCCESS, message);
}

void Logger::success(const String message)
{
    log(SUCCESS, message);
}

void Logger::warning(const char *message)
{
    log(WARNING, message);
}

void Logger::warning(const String message)
{
    log(WARNING, message);
}

void Logger::error(const char *message)
{
    log(ERROR, message);
}

void Logger::error(const String message)
{
    log(ERROR, message);
}

void Logger::critical(const char *message)
{
    log(CRITICAL, message);
}

void Logger::critical(const String message)
{
    log(CRITICAL, message);
}

void Logger::trace(const char *module, const char *message)
{
    log(TRACE, module, message);
}

void Logger::trace(const String module, const String message)
{
    log(TRACE, module, message);
}

void Logger::debug(const char *module, const char *message)
{
    log(DEBUG, module, message);
}

void Logger::debug(const String module, const String message)
{
    log(DEBUG, module, message);
}

void Logger::info(const char *module, const char *message)
{
    log(INFO, module, message);
}

void Logger::info(const String module, const String message)
{
    log(INFO, module, message);
}

void Logger::success(const char *module, const char *message)
{
    log(SUCCESS, module, message);
}

void Logger::success(const String module, const String message)
{
    log(SUCCESS, module, message);
}

void Logger::warning(const char *module, const char *message)
{
    log(WARNING, module, message);
}

void Logger::warning(const String module, const String message)
{
    log(WARNING, module, message);
}

void Logger::error(const char *module, const char *message)
{
    log(ERROR, module, message);
}

void Logger::error(const String module, const String message)
{
    log(ERROR, module, message);
}

void Logger::critical(const char *module, const char *message)
{
    log(CRITICAL, module, message);
}

void Logger::critical(const String module, const String message)
{
    log(CRITICAL, module, message);
}

void Logger::log(Level level, const char *message)
{
    log(level, "", message);
}

void Logger::log(Level level, const String message)
{
    log(level, "", message);
}

void Logger::log(Level level, const char *module, const char *message)
{
    if (level >= getLogLevel())
    {
        defaultLog(level, module, message);
    }
}

void Logger::log(Level level, const String module, const String message)
{
    log(level, module.c_str(), message.c_str());
}

String Logger::levelString(Level level)
{
    return LOG_LEVEL_STRINGS[level];
}

void Logger::defaultLog(Level level, const char *module, const char *message)
{
    if (_serial)
    {
        _serial->print(F("["));
        _serial->print(levelString(level));
        _serial->print(F("] "));

        if (strlen(module) > 0)
        {
            _serial->print(F(": "));
            _serial->print(module);
            _serial->print(F(" "));
        }

        _serial->println(message);
    }
}
