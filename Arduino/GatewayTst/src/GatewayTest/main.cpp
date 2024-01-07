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

// create the node instance on the EU-868 band
// using the radio module and the encryption key
// make sure you are using the correct band
// based on your geographical location!
LoRaWANNode node(&radio, &EU433);

float control_freq = 433.5;
float exp_freq = 433.5;
float bw = 125.0;         // Bandwidth
uint8_t control_sf = 7;   // Spreading factor for controlling
uint8_t exp_sf = 7;      // Spreading factor of experiment
uint8_t cr = 5;           // Coding rate
uint8_t syncword = 0x34;    // Syncword
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
    int state = radio.begin(control_freq, bw, exp_sf, cr, syncword, power, prelen);
    // check if module started successfully
    if (state == RADIOLIB_ERR_NONE) {
        Serial.println(F("success!"));
    } else {
        Serial.print(F("failed, code "));
        Serial.println(state);
        while (true);
    }

    // device address - this number can be anything
    // when adding new end device in TTN, you can generate this number,
    // or you can set any value you want, provided it is unique
    uint32_t devAddr = 0x12345678;

    // select some encryption keys which will be used to secure the communication
    // there are two of them - network key and application key
    // because LoRaWAN uses AES-128, the key MUST be 16 bytes (or characters) long

    // network key is the ASCII string "topSecretKey1234"
    uint8_t nwkSKey[] = { 0x74, 0x6F, 0x70, 0x53, 0x65, 0x63, 0x72, 0x65,
                            0x74, 0x4B, 0x65, 0x79, 0x31, 0x32, 0x33, 0x34 };

    // application key is the ASCII string "aDifferentKeyABC"
    uint8_t appSKey[] = { 0x61, 0x44, 0x69, 0x66, 0x66, 0x65, 0x72, 0x65,
                            0x6E, 0x74, 0x4B, 0x65, 0x79, 0x41, 0x42, 0x43 };

    // start the device by directly providing the encryption keys and device address
    Serial.print(F("[LoRaWAN] Attempting over-the-air activation ... "));
    state = node.beginABP(devAddr, nwkSKey, appSKey);
    if(state == RADIOLIB_ERR_NONE) {
        Serial.println(F("success!"));
    } else {
        Serial.print(F("failed, code "));
        Serial.println(state);
        while(true);
    }

    // // set implicit mode, use for receiving control message
    // if (radio.implicitHeader(17) == RADIOLIB_ERR_NONE) {
    //     Serial.println(F("[SX1278] 'receive' set implicitHeader, success!"));
    // } else {
    //     Serial.print(F("[SX1278] 'receive' set implicitHeader, failed"));
    //     while (true);
    // }
    // // disable CRC
    // if (radio.setCRC(false, false) == RADIOLIB_ERR_INVALID_CRC_CONFIGURATION) {
    //     Serial.println(F("[SX1278] disable CRC, failed!"));
    //     while (true);
    // } else {
    //     Serial.println(F("[SX1278] disable CRC, success!"));
    // }

}

// counter to keep track of transmitted packets
int count = 0;

void loop() {
    // send uplink to port 10
    Serial.print(F("[LoRaWAN] Sending uplink packet ... "));
    String strUp = "Hello World! #" + String(count++);
    String strDown;
    int state = node.sendReceive(strUp, 10, strDown);
    if(state == RADIOLIB_ERR_NONE) {
        Serial.println(F("received a downlink!"));

        // print data of the packet (if there are any)
        Serial.print(F("[LoRaWAN] Data:\t\t"));
        if(strDown.length() > 0) {
        Serial.println(strDown);
        } else {
        Serial.println(F("<MAC commands only>"));
        }

        // print RSSI (Received Signal Strength Indicator)
        Serial.print(F("[LoRaWAN] RSSI:\t\t"));
        Serial.print(radio.getRSSI());
        Serial.println(F(" dBm"));

        // print SNR (Signal-to-Noise Ratio)
        Serial.print(F("[LoRaWAN] SNR:\t\t"));
        Serial.print(radio.getSNR());
        Serial.println(F(" dB"));

        // print frequency error
        Serial.print(F("[LoRaWAN] Frequency error:\t"));
        Serial.print(radio.getFrequencyError());
        Serial.println(F(" Hz"));

    } else if(state == RADIOLIB_ERR_RX_TIMEOUT) {
        Serial.println(F("no downlink!"));

    } else {
        Serial.print(F("failed, code "));
        Serial.println(state);
    }

    // wait before sending another packet
    delay(30000);
}

//  ButtonPressed() is called whenever a falling edge
void ButtonPressed()
{
    Send_Flag = !Send_Flag;
}
