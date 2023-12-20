/*----------------------------------------------------------------------*
 * This Library is modified based on M5Stack I2C Common Library v1.0    *
 *                                                                      *
 * Origin work is licensed under the GNU Lesser General Public          *
 * License v2.1                                                         *
 * https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html           *
 *----------------------------------------------------------------------*/
#include "WireUtil.h"

WireUtil::WireUtil(uint8_t num)
{
    if (num == 0)
    {
        _wire = &Wire;
    }
    else
    {
        _wire = &Wire1;
    }
}

bool WireUtil::begin(int sda, int scl, uint32_t frequency)
{
    return (_wire->begin(sda, scl, frequency));
}

bool WireUtil::begin(uint8_t slaveAddr, int sda, int scl, uint32_t frequency)
{
    return (_wire->begin(slaveAddr, sda, scl, frequency));
}

bool WireUtil::writeCommand(uint8_t address, uint8_t subAddress)
{
    bool function_result = false;
    Wire.beginTransmission(address);                 // Initialize the Tx buffer
    Wire.write(subAddress);                          // Put slave register address in Tx buffer
    function_result = (Wire.endTransmission() == 0); // Send the Tx buffer

    return (function_result);
}

bool WireUtil::writeByte(uint8_t address, uint8_t subAddress, uint8_t data)
{
    bool function_result = false;
    Wire.beginTransmission(address);                 // Initialize the Tx buffer
    Wire.write(subAddress);                          // Put slave register address in Tx buffer
    Wire.write(data);                                // Put data in Tx buffer
    function_result = (Wire.endTransmission() == 0); // Send the Tx buffer

    return (function_result);
}

// Wire.h read and write protocols
bool WireUtil::writeBytes(uint8_t address, uint8_t subAddress, uint8_t *data, uint8_t length)
{
    bool function_result = false;

    Wire.beginTransmission(address);
    Wire.write(subAddress);
    for (int i = 0; i < length; i++)
    {
        Wire.write(*(data + i));
    }
    function_result = (Wire.endTransmission() == 0);

    return function_result;
}

bool WireUtil::readByte(uint8_t address, uint8_t *result)
{
    if (Wire.requestFrom(address, (uint8_t)1))
    {
        *result = Wire.read();
        return true;
    }
    return false;
}

bool WireUtil::readByte(uint8_t address, uint8_t subAddress, uint8_t *result)
{
    Wire.beginTransmission(address); // Initialize the Tx buffer
    Wire.write(subAddress);          // Put slave register address in Tx buffer
    if (Wire.endTransmission(false) == 0 &&
        Wire.requestFrom(address, (uint8_t)1))
    {
        *result = Wire.read(); // Fill Rx buffer with result
        return true;
    }
    return false;
}

bool WireUtil::readBytes(uint8_t address, uint8_t subAddress, uint8_t count, uint8_t *dest)
{
    Wire.beginTransmission(address); // Initialize the Tx buffer
    Wire.write(subAddress);          // Put slave register address in Tx buffer
    uint8_t i = 0;
    if (Wire.endTransmission(false) == 0 &&
        Wire.requestFrom(address, (uint8_t)count))
    {
        while (Wire.available())
        {
            dest[i++] = Wire.read(); // Put read results in the Rx buffer
        }
        return true;
    }
    return false;
}

bool WireUtil::readBytes(uint8_t address, uint8_t count, uint8_t *dest)
{
    uint8_t i = 0;
    if (Wire.requestFrom(address, (uint8_t)count))
    {
        while (Wire.available())
        {
            dest[i++] = Wire.read();
        }
        return true;
    }
    return false;
}

void WireUtil::scanID(bool *result)
{
    for (int i = 0x00; i <= 0x7f; i++)
    {
        *(result + i) = writeCommand(i, 0x00);
    }
}
