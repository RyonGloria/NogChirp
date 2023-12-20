#ifdef ARDUINO_USB_CDC_ON_BOOT
#undef ARDUINO_USB_CDC_ON_BOOT
#endif

// include the library
#include <Arduino.h>
#include <RadioLib.h>
#include <FastLED.h>
#include <U8g2lib.h>
#include <Wire.h>
#include <SPI.h>
#include "utils.h"

// Array of leds
CRGB leds[NUM_LEDS];

// OLED
U8G2_SSD1306_128X64_NONAME_F_HW_I2C u8g2(U8G2_R2);

// Button
void ButtonPressed();

/**
 * @brief SX1278 has the following connections:
 * @param NSS pin:   10
 * @param DIO0 pin:  2
 * @param RESET pin: 9
 * @param DIO1 pin:  3
 */
SX1278 radio = new Module(10, 15, -1, 16);

float control_freq = 433.5;   // Channel frequency for controlling
float exp_freq = 435.5;   // Channel 2 frequency
float bw = 125.0;         // Bandwidth
uint8_t control_sf = 7;   // Spreading factor for controlling
uint8_t exp_sf = 10;      // Spreading factor of experiment
uint8_t cr = 5;           // Coding rate
uint8_t syncword = 18;    // Syncword
int8_t power = 10;        // Power
uint8_t prelen = 8;       // Preamble length
float symbol_time = pow(2, exp_sf - 7); // uint ms
float esti_symbol = 30;   // estimated symbol number
int delay_max = esti_symbol * symbol_time;
boolean Send_Flag = false; // Send Flag
String data = "helloworldloraexp";  // Bin value: [810, 1010, 386, 614, 850, 406, 126, 1022, 367, 677, 347, 284, 27, 551, 991, 990, 49, 827, 822, 85, 265, 587, 421]

void setup() {
    // initialize serial port
    Serial.begin(115200);
    delay(200);

    pinMode(USER_BTN, INPUT);   // set USER_BTN as input
    attachInterrupt(USER_BTN, ButtonPressed, FALLING); // attach interrupt to USER_BTN, FALLING mode

    // reset the sx1278 radio through IO expander
    Wire.begin();
    IOreset();

    // initialize OLED
    u8g2.begin();
    u8g2.clearBuffer();
    delay(200);

    // initialize LED
    FastLED.addLeds<NEOPIXEL, RGB_BUILTIN>(leds, NUM_LEDS); // GRB ordering is assumed
    FastLED.setBrightness(RGB_BRIGHTNESS / 8);              // Set brightness (0-255)
    FastLED.clear(true);                                    // clear all pixel data
    delay(200);

    // initialize SX1278 with default settings
    Serial.print(F("[SX1278] Initializing ... "));
    // @params frequency, bandwidth, spreading factor, coding rate, syncword, power, preamble length
    int state = radio.begin(control_freq, bw, control_sf, cr, syncword, power, prelen);
    // check if module started successfully
    if (state == RADIOLIB_ERR_NONE) {
        Serial.println(F("success!"));
    } else {
        Serial.print(F("failed, code "));
        Serial.println(state);
        while (true);
    }

    // set implicit mode, use for receiving control message
    if (radio.implicitHeader(17) == RADIOLIB_ERR_NONE) {
        Serial.println(F("[SX1278] 'receive' set implicitHeader, success!"));
    } else {
        Serial.print(F("[SX1278] 'receive' set implicitHeader, failed"));
        while (true);
    }

    // disable CRC
    if (radio.setCRC(false, false) == RADIOLIB_ERR_INVALID_CRC_CONFIGURATION) {
        Serial.println(F("[SX1278] disable CRC, failed!"));
        while (true);
    } else {
        Serial.println(F("[SX1278] disable CRC, success!"));
    }


}

void loop() {
    // wait for incoming control message
    String str;
    if (radio.receive(str) == RADIOLIB_ERR_NONE) {
        // print the data of the packet
        Serial.print(F("[SX1278] Data:\t"));
        Serial.println(str);

        // judge if the control message is "Start"
        if (str.indexOf("Start") != -1) {
            Serial.println(F("[SX1278] get the Right command!"));

            // set carrier frequency to exp_freq
            if (radio.setFrequency(exp_freq) == RADIOLIB_ERR_INVALID_FREQUENCY) {
                Serial.println(F("[SX1278] set to exp_freq, failed!"));
                while (true);
            } else {
                Serial.println(F("[SX1278] set to exp_freq, success!"));
            }

            // set spreading factor to exp_sf
            if (radio.setSpreadingFactor(exp_sf) == RADIOLIB_ERR_INVALID_SPREADING_FACTOR) {
                Serial.println(F("[SX1278] set to exp_sf, failed!"));
                while (true);
            } else {
                Serial.println(F("[SX1278] set to exp_sf, success!"));
            }

            // set implicit mode, use for sending experiment message
            if (radio.implicitHeader(sizeof(data)) == RADIOLIB_ERR_NONE) {
                Serial.println(F("[SX1278] 'transmit' set implicitHeader, success!"));
            } else {
                Serial.print(F("[SX1278] 'transmit' set implicitHeader, failed"));
                while (true);
            }

            Send_Flag = true; // set Send_Flag to true

        } else {
            Serial.println(F("[SX1278] get the Wrong command!"));
        }
    }

    // if receive the control message, start sending experiment message
    if (Send_Flag)
    {
        // transmit packet
        int state = radio.transmit(data);
        if (state == RADIOLIB_ERR_NONE) {
            Serial.println("transmit with " + String(exp_freq) + "MHz, success!");
        }

        leds[0] = CRGB::DarkRed;
        FastLED.show();
        delay(1000);
        FastLED.showColor(CRGB(0,0,0), 0);
        delay(1000);
    }
}

//  ButtonPressed() is called whenever a falling edge
void ButtonPressed()
{
    Send_Flag = !Send_Flag;
}
