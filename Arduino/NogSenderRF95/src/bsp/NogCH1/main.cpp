/*
 * LEMoN BSP Demo
 *
 * Project Config:
 *     src_dir = src/bsp
 *
 * Use env: default
 */

#include <Arduino.h>

#include <Lemon.h>
#include <Logger.h>
#include <Button.h>
#include <Oled.h>
#include <Led.h>

#include <RH_RF95.h>

void ButtonPressed();
void oled_display(String msg, String msg2 = "", String msg3 = "", String msg4 = "");

Lemon lemon;

// Button btn = Button(USER_BTN, true, 10);
Logger logger;
Oled oled;
Led led;
boolean Send_Flag = false;

RH_RF95 rf95(RA_CS, DIO0);
float frequency = 435.96875;  // 434.90625MHz ~ 435.03125MHz

void setup()
{
    // put your setup code here, to run once:

    // Logger init
    // Output at USB2
    logger.setLogLevel(INFO);
    logger.begin(115200);
    delay(200);

    pinMode(USER_BTN, INPUT);   //按键设置为输入模式
    attachInterrupt(USER_BTN, ButtonPressed, FALLING); //按键中断设置为下降沿触发

    // Board init
    // This will reset LoRa and Camera
    lemon.init();
    logger.info("LEMoN begin!");

    // RF95 init
    if (!rf95.init())
    {
        logger.error("RF95 init failed!");
        while (1)
            ;
    }

    // LoRa parameters Setup
    rf95.setFrequency(437.0);                 // Setup ISM frequency
    rf95.setTxPower(13);                    // Setup Power,dBm
    rf95.setSpreadingFactor(10);             // Setup Spreading Factor (6 ~ 12)
    rf95.setSignalBandwidth(125000);        // Setup BandWidth, option: 7800,10400,15600,20800,31200,41700,62500,125000,250000,500000（Lower BandWidth for longer distance）
    rf95.setCodingRate4(5);                 // Setup Coding Rate:5(4/5),6(4/6),7(4/7),8(4/8)

    // OLED display
    oled.init();
    oled.setCursor(0, 15);
    oled.print("Initializing...");
    oled.sendBuffer();

    // RGB LED display
    led.init();
    led.showColor(CRGB::Red);
    delay(500);
    led.showColor(CRGB::Green);
    delay(500);
    led.showColor(CRGB::White);
    delay(500);

    // Init finished print
    logger.success("Init finished.");
    oled.clearBuffer();
    oled.setCursor(0, 15);
    oled.print("Init finished.");
    oled.sendBuffer();
}

//  ButtonPressed()中断函数
void ButtonPressed()
{
    Send_Flag = !Send_Flag;
}

void loop()
{
    // put your main code here, to run repeatedly:

    // Wait for Msg from Controller
    if (rf95.available())  // Will be set by the interrupt handler when a good message is received
    {
        // Should be a message for us now
        uint8_t buf[RH_RF95_MAX_MESSAGE_LEN];
        uint8_t len = sizeof(buf);
        if (rf95.recv(buf, &len))
        {
            logger.info("got reply: ");
            logger.info((char*)buf);
            // oled_display((char *)buf, 0, 15);
            Send_Flag = !Send_Flag;
            rf95.setFrequency(frequency);  // Setup frequency of 435.96875MHz
        }else{
            logger.info("Wait for Msg from Controller.");
        }
    }

    // Send a message, Device 1 use frequency of 435.96875MHz
    if (Send_Flag)
    {
        uint8_t data[] = "helloworld";
        rf95.setPayloadCRC(false); // Disable CRC
        rf95.sendWithoutHeader(data, sizeof(data));
        rf95.waitPacketSent();
        oled_display("Send with " + String(frequency) + "MHz");
        logger.info("Sent a msg : 'A 20-byte payload...'");
        led.showColor(CRGB::Green);
        // Send_Flag = !Send_Flag;
    }else{
        oled_display("Waiting...");
        led.showColor(CRGB::Red);
    }

    delay(1000);
}


/***
 * @brief OLED显示函数
 * @author by ChenKai
 * @param msg: 第一行显示的字符串 @param msg2: 第二行显示的字符串 @param msg3: 第三行显示的字符串 @param msg4: 第四行显示的字符串
 * @note 如果字符串长度超过了 20 个字符，会自动换行显示
*/
void oled_display(String msg, String msg2, String msg3, String msg4)
{
    const uint8_t MAX_LINES = 4, MAX_CHARS_PER_LINE = 20;  // 最大行数和每行最大字符数
    String messages[MAX_LINES] = {msg, msg2, msg3, msg4};  // 字符串数组
    uint8_t line_num[MAX_LINES] = {15, 30, 45, 60};        // 每行的 y 坐标
    bool flag = false;
    uint8_t index = 0, begin = 0;  // begin 用于记录起始字符

    oled.clearBuffer(); // 清除缓存

    for (uint8_t line = 0; line < MAX_LINES; line++)  // 遍历每一行，每一行限定 20 个字符
    {
        uint8_t start = begin++ * MAX_CHARS_PER_LINE;   // 每一行的起始字符
        uint8_t end = start + MAX_CHARS_PER_LINE;
        if (end >= messages[index].length())   // 如果超过了字符串长度
        {
            end = messages[index].length(); flag = true;
        }

        oled.setCursor(0, line_num[line]);
        oled.print(messages[index].substring(start, end));
        if (flag)  // 如果超过了字符串长度，切换到下一个需要换行显示的字符串
        {
            index++; begin = 0; flag = false;
        }

    }
    oled.sendBuffer(); // 发送缓存
}
