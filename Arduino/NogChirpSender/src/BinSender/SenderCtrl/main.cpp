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
#include <cmath>

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

// flag to indicate that a packet was received
volatile bool transmittedFlag = false;

// flag to indicate frequency must be changed
volatile bool fhssChangeFlag = false;

float control_freq = 433.5;   // Channel frequency for controlling
float exp_freq = 433.5;       // Channel 2 frequency
float bw = 125.0;             // Bandwidth
uint8_t control_sf = 7;       // Spreading factor for controlling
uint8_t exp_sf = 10;          // Spreading factor of experiment
uint8_t cr = 5;               // Coding rate
uint8_t syncword = 18;        // Syncword
int8_t power = 10;            // Power
uint8_t prelen = 8;           // Preamble length
float symbol_time = pow(2, exp_sf - 7);    // uint ms
float esti_symbol = 30;       // estimated symbol number
int delay_max = esti_symbol * symbol_time;
boolean Send_Flag = false;    // Send Flag
int true_bin[23] = {0};
int true_bin_SF10[23] = {809, 1009, 385, 613, 849, 405, 125, 1021, 366, 676, 346, 283, 26, 550, 990, 989, 48, 826, 53, 493, 264, 583, 420};
int bin = 2<<(exp_sf-1);      // bin range = 2^sf


int transmissionState = RADIOLIB_ERR_NONE;   // save transmission state between loops
int true_bin_count = 0;                      // record the number of successful bins

// 对 LoRa 解码的格雷和汉明进行反编码得到可以使 chirpbin 固定的发送序列
uint8_t data[] = {0xFF, 0xFE, 0xFC, 0xF8, 0xF0, 0xE1, 0xC2, 0x85, 0x0B, 0x17, 0x2F, 0x5E, 0xBC, 0x78, 0xF1, 0xE3,
                  0xC6, 0x8D, 0x1A, 0x34, 0x68, 0xD0, 0xA0, 0x40, 0x80, 0x01, 0x02, 0x04, 0x08, 0x11, 0x23, 0x47,
                  0x8E, 0x1C, 0x38, 0x71, 0xE2, 0xC4, 0x89, 0x12, 0x25, 0x4B, 0x97, 0x2E, 0x5C, 0xB8, 0x70, 0xE0,
                  0xC0, 0x81, 0x03, 0x06, 0x0C, 0x19, 0x32, 0x64, 0xC9, 0x92, 0x24, 0x49, 0x93, 0x26, 0x4D, 0x9B,
                  0x37, 0x6E, 0xDC, 0xB9, 0x72, 0xE4, 0xC8, 0x90, 0x20, 0x41, 0x82, 0x05, 0x0A, 0x15, 0x2B, 0x56};
int dataLength;

// this function is called when a complete packet
// is transmitted by the module
// IMPORTANT: this function MUST be 'void' type and MUST NOT have any arguments!
#if defined(ESP8266) || defined(ESP32)
  ICACHE_RAM_ATTR
#endif
void setTxFlag(void) {
  transmittedFlag = true;
}

// this function is called when FhssChangeChannel interrupt occurs
// (at the beginning of each transmission)
// IMPORTANT: this function MUST be 'void' type and MUST NOT have any arguments!
#if defined(ESP8266) || defined(ESP32)
  ICACHE_RAM_ATTR
#endif
void setFHSSFlag(void) {
  fhssChangeFlag = true;
}

void setup() {
    dataLength = 23;
    // set true_bin
    if(exp_sf == 10){
        for(int i = 0; i < 23; i++)
            true_bin[i] = true_bin_SF10[i] - 1;
    }

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
    Serial.print(F("[SX1278] Initializing... "));
    // @params frequency, bandwidth, spreading factor, coding rate, syncword, power, preamble length
    int state = radio.begin(control_freq, bw, control_sf, cr, syncword, power, prelen);
    // check if module started successfully
    if (state == RADIOLIB_ERR_NONE) {
        Serial.println(F("success!"));
    } else {
        Serial.print(F(", failed, code "));
        Serial.println(state);
        while (true);
    }
    delay(200);

    // set hop period in symbols, this will also enable FHSS
    Serial.print(F("[SX1278] set FHSSHoppingPeriod in symbols... "));
    state = radio.setFHSSHoppingPeriod(1);
    if (state == RADIOLIB_ERR_NONE) {
        Serial.println(F("success!"));
    } else {
        Serial.print(F("failed, code "));
        Serial.println(state);
        while (true);
    }
    delay(200);

    // set the function to call when transmission is finished
    radio.setDio0Action(setTxFlag, RISING);
    // set the function to call when we need to change frequency
    radio.setDio1Action(setFHSSFlag, RISING);

    // set implicit mode
    Serial.print(F("[SX1278] set implicitHeader... "));
    state = radio.implicitHeader(dataLength);
    if (state == RADIOLIB_ERR_NONE) {
        Serial.println(F("success!"));
    } else {
        Serial.print(F("failed, code "));
        Serial.println(state);
        while (true);
    }

    // disable CRC
    Serial.print(F("[SX1278] disable CRC... "));
    state = radio.setCRC(false, false);
    if (state == RADIOLIB_ERR_INVALID_CRC_CONFIGURATION) {
        Serial.print(F("failed, code "));
        Serial.println(state);
        while (true);
    } else {
        Serial.println(F("success!"));
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
        // detect if the interrupt is triggered
        if (transmittedFlag == true) {
          // reset flag
          transmittedFlag = false;
          // check if the transmission was successful
          if (transmissionState == RADIOLIB_ERR_NONE) {
            Serial.println(F("finished!"));
          } else {
            Serial.print(F("failed, code "));
            Serial.println(transmissionState);
          }

          // reset parameters
          true_bin_count = 0;

          // return to home channel before the next transaction
          radio.setFrequency(exp_freq);

          // wait seconds before transmitting again
          delay(6000);
          leds[0] = CRGB::DarkRed;
          FastLED.show();
          delay(1000);
          FastLED.showColor(CRGB(0,0,0), 0);

          // 发送下一个数据包
          Serial.println(F("[SX1278] Start transferring packet..."));
          Serial.print(F("         "));
          transmissionState = radio.startTransmit(data, dataLength, 0U);
        }

        // deal with FHSS interrupt
        if (fhssChangeFlag == true) {
          // change the frequency according to the true_bin
          float adjust_bw = bw/1000 * (float)true_bin[true_bin_count]/bin;
          // set carrier frequency
          // if (radio.setFrequency(exp_freq + adjust_bw) == RADIOLIB_ERR_INVALID_FREQUENCY) {
          //     Serial.print("[SX1278] freq set to " + String(exp_freq + adjust_bw) + ", Bin Value (" + String(true_bin[true_bin_count]) +"), failed!");
          // } else {
          //    Serial.print("[SX1278] freq set to " + String(exp_freq + adjust_bw) + ", Bin Value (" + String(true_bin[true_bin_count]) +"), success!");
          // }
          if (radio.setFrequency(exp_freq + adjust_bw) != RADIOLIB_ERR_INVALID_FREQUENCY) {
              Serial.print(String(true_bin[true_bin_count]) + " ");
          }
          true_bin_count++;

          // clear the FHSS interrupt
          radio.clearFHSSInt();

          // we're ready to do another hop, clear the flag
          fhssChangeFlag = false;
        }
    }

}

//  ButtonPressed() is called whenever a falling edge
void ButtonPressed()
{
    Send_Flag = !Send_Flag;
}
