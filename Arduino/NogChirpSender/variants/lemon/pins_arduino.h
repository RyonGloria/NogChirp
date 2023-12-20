#ifndef Pins_Arduino_h
#define Pins_Arduino_h

#include <stdint.h>

#define USB_VID 0x303a
#define USB_PID 0x1001

#define EXTERNAL_NUM_INTERRUPTS 46
#define NUM_DIGITAL_PINS 48
#define NUM_ANALOG_INPUTS 18

static const uint8_t LED_BUILTIN = 45;
#define BUILTIN_LED LED_BUILTIN
#define LED_BUILTIN LED_BUILTIN
#define RGB_BUILTIN LED_BUILTIN
#define NUM_LEDS 1

#define USER_BTN 0

#define SD_CS 9

#define RA_CS 10
#define DIO0 15
#define DIO1 16

#define I2C_SDA 17
#define I2C_SCL 18

#define OV2640_ADDR 0x30   // Camera
#define SSD1306_ADDR 0x34  // OLED Screen
#define PCA9554A_ADDR 0x38 // IO Extender

#define PIN_PWDN -1
#define PIN_RESET -1
#define PIN_XCLK 40
#define PIN_SIOD 17
#define PIN_SIOC 18
#define PIN_D7 39
#define PIN_D6 41
#define PIN_D5 42
#define PIN_D4 12
#define PIN_D3 3
#define PIN_D2 14
#define PIN_D1 47
#define PIN_D0 13
#define PIN_VSYNC 21
#define PIN_HREF 38
#define PIN_PCLK 11

#define analogInputToDigitalPin(p) (((p) < 20) ? (analogChannelToDigitalPin(p)) : -1)
#define digitalPinToInterrupt(p) (((p) < 48) ? (p) : -1)
#define digitalPinHasPWM(p) (p < 46)

static const uint8_t TX = 43;
static const uint8_t RX = 44;

static const uint8_t SDA = 17;
static const uint8_t SCL = 18;

static const uint8_t SS = 10;
static const uint8_t MOSI = 5;
static const uint8_t MISO = 4;
static const uint8_t SCK = 48;

static const uint8_t A0 = 1;
static const uint8_t A1 = 2;
static const uint8_t A2 = 3;
static const uint8_t A3 = 4;
static const uint8_t A4 = 5;
static const uint8_t A5 = 6;
static const uint8_t A6 = 7;
static const uint8_t A7 = 8;
static const uint8_t A8 = 9;
static const uint8_t A9 = 10;
static const uint8_t A10 = 11;
static const uint8_t A11 = 12;
static const uint8_t A12 = 13;
static const uint8_t A13 = 14;
static const uint8_t A14 = 15;
static const uint8_t A15 = 16;
static const uint8_t A16 = 17;
static const uint8_t A17 = 18;

static const uint8_t T1 = 1;
static const uint8_t T2 = 2;
static const uint8_t T3 = 3;
static const uint8_t T4 = 4;
static const uint8_t T5 = 5;
static const uint8_t T6 = 6;
static const uint8_t T7 = 7;
static const uint8_t T8 = 8;
static const uint8_t T9 = 9;
static const uint8_t T10 = 10;
static const uint8_t T11 = 11;
static const uint8_t T12 = 12;
static const uint8_t T13 = 13;
static const uint8_t T14 = 14;

#endif /* Pins_Arduino_h */
