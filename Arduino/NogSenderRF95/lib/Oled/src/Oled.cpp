#include "Oled.h"

Oled::Oled() : U8G2_SSD1306_128X64_NONAME_F_HW_I2C(U8G2_R2)
{
}

void Oled::init()
{
    begin();
    enableUTF8Print();
    setFont(u8g2_font_wqy12_t_chinese3);
    setFontDirection(0);
    clearBuffer();
}

Oled OLED;
