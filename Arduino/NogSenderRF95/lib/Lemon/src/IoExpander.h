#ifndef _IOEXPANDER_H_
#define _IOEXPANDER_H_

#include <Arduino.h>
#include "WireUtil.h"

// PCA9554 Command Byte
#define INPUTPORT 0x00
#define OUTPUTPORT 0x01
#define POLINVPORT 0x02
#define CONFIGPORT 0x03

#define ALLOUTPUT 0x00
#define ALLINPUT 0xFF

class IoExpander
{
public:
    IoExpander(uint8_t SlaveAddress);

    bool twiRead(uint8_t address, uint8_t *data);
    bool twiWrite(uint8_t address, uint8_t data);

    bool pinMode(uint8_t pin, bool state);
    bool portMode(uint8_t value);

    bool digitalWrite(uint8_t pin, bool state);
    bool digitalWritePort(uint8_t value);

    // Return 0xFF when failed
    uint8_t digitalRead(uint8_t pin);
    uint8_t digitalReadPort();

    // Return status to status parameter
    uint8_t digitalReadWithStatus(uint8_t pin, bool *status);
    uint8_t digitalReadPortWithStatus(bool *status);

private:
    int _addr;
    WireUtil _wire = WireUtil(0);
};

extern IoExpander PIO;

#endif
