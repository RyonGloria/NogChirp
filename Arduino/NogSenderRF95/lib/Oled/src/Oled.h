#ifndef _OLED_H_
#define _OLED_H_

#include <Arduino.h>
#include <U8g2lib.h>
#include <Wire.h>

class Oled : public U8G2_SSD1306_128X64_NONAME_F_HW_I2C
{
public:
    Oled();
    void init();

private:
};

extern Oled OLED;

#warning "OLED library is WIP."

#endif
