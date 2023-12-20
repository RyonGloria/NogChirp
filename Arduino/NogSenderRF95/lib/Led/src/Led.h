#ifndef _LED_H_
#define _LED_H_

#include <Arduino.h>
#include <FastLED.h>

class Led
{
public:
    Led();

    void init();
    void setBrightness(uint8_t scale);
    void showColor(CRGB color);
    void showColor(uint32_t color);

private:
};

#endif
