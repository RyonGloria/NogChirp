/**
 * \par Copyright (C), 2022-2023, Mobinets
 * \brief   LEMoN library.
 * @file    Lemon.h
 * @version V0.0.1
 * @date    2023/04/20
 * @brief   Header for Lemon.cpp module
 *
 * \par Description
 * This file is a drive for LEMoN core.
 *
 * \par Method List:
 *
 *  System:
        lemon.begin();
 */

#ifndef _LEMON_H_
#define _LEMON_H_

#if defined(ESP32)

#include <Arduino.h>

class Lemon
{
public:
    Lemon();
    void init();
    void resetLoRa();
    void resetCamera();

private:
    bool isInited;
};

#else
#error "This library only supports boards with ESP32 processor."
#endif
#endif
