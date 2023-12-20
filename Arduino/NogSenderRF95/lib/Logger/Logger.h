#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <Arduino.h>

enum Level
{
    TRACE = 0,
    DEBUG,
    INFO,
    SUCCESS,
    WARNING,
    ERROR,
    CRITICAL
};

class Logger
{
public:
    Logger();

    void begin(uint32_t baud = 115200);

    void setSerial(HardwareSerial &serial);
    HardwareSerial &getSerial();

    void setLogLevel(Level level);
    Level getLogLevel();

    void trace(const char *message);
    void trace(const String message);
    void debug(const char *message);
    void debug(const String message);
    void info(const char *message);
    void info(const String message);
    void success(const char *message);
    void success(const String message);
    void warning(const char *message);
    void warning(const String message);
    void error(const char *message);
    void error(const String message);
    void critical(const char *message);
    void critical(const String message);

    void trace(const char *module, const char *message);
    void trace(const String module, const String message);
    void debug(const char *module, const char *message);
    void debug(const String module, const String message);
    void info(const char *module, const char *message);
    void info(const String module, const String message);
    void success(const char *module, const char *message);
    void success(const String module, const String message);
    void warning(const char *module, const char *message);
    void warning(const String module, const String message);
    void error(const char *module, const char *message);
    void error(const String module, const String message);
    void critical(const char *module, const char *message);
    void critical(const String module, const String message);

    void log(Level level, const char *message);
    void log(Level level, const String message);

    void log(Level level, const char *module, const char *message);
    void log(Level level, const String module, const String message);

    String levelString(Level level);

private:
    void defaultLog(Level level, const char *module, const char *message);

    Level _level = INFO;
    HardwareSerial *_serial = &Serial;
};

#endif // _LOGGER_H_
