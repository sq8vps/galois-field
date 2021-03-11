#ifndef GF_H
#define GF_H

#include <stdint.h>
#include <string.h>



class GF
{
public:
    uint16_t add(uint16_t x, uint16_t y);
    uint16_t sub(uint16_t x, uint16_t y);
    uint16_t mul(uint16_t x, uint16_t y);
    uint16_t div(uint16_t dividend, uint16_t divisor);
    uint16_t pow(uint16_t x, uint16_t exponent);
    uint16_t inv(uint16_t x);
    uint16_t slowMul(uint16_t x, uint16_t y);

    GF(uint16_t p);
    ~GF();



private:
    uint16_t *exp;
    uint16_t *log;
    uint16_t len;
};


#endif // REEDSOLOMON_H
