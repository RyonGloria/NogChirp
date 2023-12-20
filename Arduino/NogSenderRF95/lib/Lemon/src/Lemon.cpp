#include <Wire.h>

#include "Lemon.h"
#include "IoExpander.h"

Lemon::Lemon() : isInited(0)
{
}

void Lemon::init()
{
    // Correct init once
    if (isInited)
    {
        return;
    }
    else
    {
        isInited = true;
    }

    // Required by IO Expander
    Wire.begin();

    // Set IO Expander
    PIO.portMode(ALLINPUT);
    PIO.pinMode(EX_PIN_CAM_RST, OUTPUT);
    PIO.pinMode(EX_PIN_RA_RST, OUTPUT);

    // Reset LoRa and Camera
    resetLoRa();
    resetCamera();
}

void Lemon::resetLoRa()
{
    PIO.digitalWrite(EX_PIN_RA_RST, LOW);
    delay(50);
    PIO.digitalWrite(EX_PIN_RA_RST, HIGH);
}

void Lemon::resetCamera()
{
    PIO.digitalWrite(EX_PIN_CAM_RST, LOW);
    delay(50);
    PIO.digitalWrite(EX_PIN_CAM_RST, HIGH);
}
