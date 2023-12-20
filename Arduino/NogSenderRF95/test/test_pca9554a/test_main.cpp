#include <Arduino.h>
#include <Wire.h>
#include <SPI.h>

#include <unity.h>

#include <PCA9554.h>

PCA9554 io(PCA9554A_ADDR);

void setUp(void)
{
    // set stuff up here
}

void tearDown(void)
{
    // clean stuff up here
}

void test_iic_pin_sda(void)
{
    TEST_ASSERT_EQUAL(17, SDA);
}

void test_iic_pin_scl(void)
{
    TEST_ASSERT_EQUAL(18, SCL);
}

void test_io_high(void)
{
    byte pin = 1;
    io.digitalWrite(pin, HIGH);
    io.digitalRead(pin);
    TEST_ASSERT_EQUAL(HIGH, pin);
}

void test_io_low(void)
{
    byte pin = 1;
    io.digitalWrite(pin, LOW);
    io.digitalRead(pin);
    TEST_ASSERT_EQUAL(LOW, pin);
}

uint8_t i = 0;
uint8_t max_tests = 5;

void setup()
{
    // NOTE!!! Wait for >2 secs
    // if board doesn't support software reset via Serial.DTR/RTS
    delay(2000);

    Wire.begin();
    io.portMode(ALLOUTPUT);
    io.pinMode(0, INPUT);

    UNITY_BEGIN(); // IMPORTANT LINE!
    delay(500);

    RUN_TEST(test_iic_pin_sda);
    delay(500);
    RUN_TEST(test_iic_pin_scl);

    while (i < max_tests)
    {
        RUN_TEST(test_io_high);
        delay(500);
        RUN_TEST(test_io_low);
        delay(500);
        i++;
    }

    UNITY_END(); // stop unit testing
}

void loop()
{
}
