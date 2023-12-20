#include "IoExpander.h"

IoExpander::IoExpander(uint8_t address)
{
    _addr = address;
}

bool IoExpander::twiRead(uint8_t address, uint8_t *data)
{
    bool function_result = false;
    function_result = _wire.readByte(_addr, address, data);
    return function_result;
}

bool IoExpander::twiWrite(uint8_t address, uint8_t data)
{
    bool function_result = false;
    function_result = _wire.writeByte(_addr, address, data);
    return function_result;
}

bool IoExpander::pinMode(uint8_t pin, bool state)
{
    uint8_t data = 0xFF;
    bool function_result = false;
    function_result = this->twiRead(CONFIGPORT, &data);
    if (function_result && (pin <= 7))
    {
        if (!state)
        {
            data |= (1 << pin);
            if (this->portMode(data))
            {
                return true;
            }
            return false;
        }
        else if (state)
        {
            data &= ~(1 << pin);
            if (this->portMode(data))
            {
                return true;
            }
            return false;
        }
    }
    return false;
}

bool IoExpander::portMode(uint8_t value)
{
    if (this->twiWrite(CONFIGPORT, value))
    {
        return true;
    }
    return false;
}

bool IoExpander::digitalWrite(uint8_t pin, bool state)
{
    uint8_t data = 0xFF;
    bool function_result = false;
    function_result = this->twiRead(OUTPUTPORT, &data);
    if (function_result && pin <= 7)
    {
        if (state)
        {
            data |= (1 << pin);
            if (this->digitalWritePort(data))
                return true;
            return false;
        }
        else if (!state)
        {
            data &= ~(1 << pin);
            if (this->digitalWritePort(data))
                return true;
            return false;
        }
    }
    return false;
}

bool IoExpander::digitalWritePort(uint8_t value)
{
    if (this->twiWrite(OUTPUTPORT, value))
        return true;
    return false;
}

uint8_t IoExpander::digitalRead(uint8_t pin)
{
    uint8_t data = 0xFF;
    bool function_result = false;
    function_result = this->twiRead(INPUTPORT, &data);
    if (function_result && (pin <= 7))
    {
        data &= (1 << pin);
        return data;
    }
    return 0xFF;
}

uint8_t IoExpander::digitalReadPort()
{
    uint8_t data = 0xFF;
    bool function_result = false;
    function_result = this->twiRead(INPUTPORT, &data);
    if (function_result)
    {
        return data;
    }
    return 0xFF;
}

uint8_t IoExpander::digitalReadWithStatus(uint8_t pin, bool *status)
{
    uint8_t data = 0xFF;
    *status = this->twiRead(INPUTPORT, &data);
    if ((*status) && (pin <= 7))
    {
        data &= (1 << pin);
        if (data > 0)
        {
            data = 1;
        }
        else
        {
            data = 0;
        }
        return data;
    }
    return 0xFF;
}

uint8_t IoExpander::digitalReadPortWithStatus(bool *status)
{
    uint8_t data = 0xFF;
    *status = this->twiRead(INPUTPORT, &data);
    if ((*status))
    {
        return data;
    }
    return 0xFF;
}

IoExpander PIO = IoExpander(PCA9554A_ADDR);
