#include "utils.h"
#include <Wire.h>

PCA9554 io(PCA9554A_ADDR);

void IOreset(void){
    io.pinMode(6,OUTPUT);
    io.pinMode(7,OUTPUT);
    io.digitalWrite(6,LOW);
    delay(1);
    io.digitalWrite(6,HIGH);
    delay(5);
    io.digitalWrite(7,LOW);
    delay(1);
    io.digitalWrite(7,HIGH);
    delay(5);
}