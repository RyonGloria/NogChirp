#include "Led.h"

// Array of leds
CRGB leds[NUM_LEDS];

Led::Led()
{
}

void Led::init()
{
    FastLED.addLeds<WS2812, RGB_BUILTIN, GRB>(leds, NUM_LEDS); // GBR Order
    FastLED.setBrightness(RGB_BRIGHTNESS / 8);
}

void Led::setBrightness(uint8_t scale)
{
    FastLED.setBrightness(scale);
}

void Led::showColor(CRGB color)
{
    leds[0] = color;
    FastLED.show();
}

void Led::showColor(uint32_t color)
{
    leds[0] = CRGB(color);
    FastLED.show();
}
