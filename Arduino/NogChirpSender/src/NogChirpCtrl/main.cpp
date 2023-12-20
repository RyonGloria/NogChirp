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

// this function is called when a complete packet
// is transmitted by the module
// IMPORTANT: this function MUST be 'void' type
//            and MUST NOT have any arguments!
#if defined(ESP8266) || defined(ESP32)
  ICACHE_RAM_ATTR
#endif

float control_freq = 433.5;       // Channel frequency for controlling
float exp_freq_ch1 = 433.46875;   // Channel 1 frequency, delay 50ms
float exp_freq_ch2 = 433.5;       // Channel 2 frequency
float exp_freq_ch3 = 433.53125;   // Channel 3 frequency, delay 100ms
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
boolean Btn_Flag = false; // Button Flag
String data = "CtrlMsg:StartSend";

void setup() {
    // initialize serial port
    Serial.begin(115200);
    delay(200);

    pinMode(USER_BTN, INPUT);   // set USER_BTN as input
    attachInterrupt(USER_BTN, ButtonPressed, FALLING); // attach interrupt to USER_BTN, FALLING mode

    Wire.begin();
    // reset the sx1278 radio through IO expander
    IOreset();
    u8g2.begin();
    u8g2.clearBuffer();
    delay(500);
    FastLED.addLeds<NEOPIXEL, RGB_BUILTIN>(leds, NUM_LEDS); // GRB ordering is assumed
    FastLED.setBrightness(RGB_BRIGHTNESS / 8);              // Set brightness (0-255)
    FastLED.clear(true);                                    // clear all pixel data
    delay(500);

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

    // set implicit mode
    if (radio.implicitHeader(sizeof(data)) == RADIOLIB_ERR_NONE) {
        Serial.println(F("[SX1278] set implicitHeader mode, success!"));
    } else {
        Serial.print(F("[SX1278] set implicitHeader mode, failed"));
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
    if (Btn_Flag)
    {
        // transmit control message
        int state = radio.transmit(data);
        if (state == RADIOLIB_ERR_NONE) {
            Serial.println("transmit '" + String(data) + "' success!");
        }

        Btn_Flag = false;
        leds[0] = CRGB::DarkGreen;
        FastLED.show();
        delay(1000);
        FastLED.showColor(CRGB(0,0,0), 0);
    }
}

//  ButtonPressed() is called whenever a falling edge
void ButtonPressed()
{
    Btn_Flag = !Btn_Flag;
}

