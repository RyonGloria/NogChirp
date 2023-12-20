/*----------------------------------------------------------------------*
 * This Library is modified based on M5Stack I2C Common Library v1.0    *
 *                                                                      *
 * Origin work is licensed under the GNU Lesser General Public          *
 * License v2.1                                                         *
 * https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html           *
 *----------------------------------------------------------------------*/
#ifndef WireUtil_h
#define WireUtil_h

#include <Arduino.h>
#include <Wire.h>

class WireUtil
{
public:
    WireUtil(uint8_t num);
    bool begin(int sda, int scl, uint32_t frequency = 0);
    bool begin(uint8_t slaveAddr, int sda, int scl, uint32_t frequency);
    inline bool begin()
    {
        return begin(-1, -1, static_cast<uint32_t>(0));
    }
    inline bool begin(uint8_t addr)
    {
        return begin(addr, -1, -1, 0);
    }
    inline bool begin(int addr)
    {
        return begin(static_cast<uint8_t>(addr), -1, -1, 0);
    }
    bool writeCommand(uint8_t address, uint8_t subAddress);
    bool writeByte(uint8_t address, uint8_t subAddress, uint8_t data);
    bool writeBytes(uint8_t address, uint8_t subAddress, uint8_t *data, uint8_t length);
    bool readByte(uint8_t address, uint8_t *result);
    bool readByte(uint8_t address, uint8_t subAddress, uint8_t *result);
    bool readBytes(uint8_t address, uint8_t count, uint8_t *dest);
    bool readBytes(uint8_t address, uint8_t subAddress, uint8_t count, uint8_t *dest);
    void scanID(bool *result);

private:
    TwoWire *_wire;
};
#endif
