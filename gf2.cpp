/*
    This file is part of simple Galois field library.

    This is a free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    It is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
* @file gf2.h
* @brief Simple Galois field library for GF(2^8)
* @version 1.1
* @author Piotr Wilkon <sq8vps@gmail.com>
* @copyright Copyright 2021 Piotr Wilkon, licensed under GNU GPLv3
**/

#include "gf2.h"


uint8_t GF2::add(uint8_t x, uint8_t y)
{
    return x ^ y; //addition is done by XOR in GF(2^8)
}
uint8_t GF2::sub(uint8_t x, uint8_t y)
{
  return x ^ y; //subtraction is done by XOR in GF(2^8) - surprisingly, addition is the same as subtraction in GF(2^n)
}
/**
 * @brief Fast multiplication in Galois field
 * @param x Multiplicand
 * @param y Multiplier
 * @return Multiplication result
 */
uint8_t GF2::mul(uint8_t x, uint8_t y)
{
    if((x == 0) || (y == 0)) //trivial multiplication by 0
        return 0;
    //fast multiplication using lookup tables
    //we know that log(x)+log(y)=log(x*y), and b^log(a)=a, when b is the logarithm base
    //so b^(log(x)+log(y))=b^log(x*y)=x*y, where b is the logarithm base
    return exp[log[x] + log[y]];
}

/**
 * @brief Fast division in Galois field
 * @param dividend Dividend
 * @param divisor Divisor
 * @return Division result. 0 is returned when dividing by 0.
 */
uint8_t GF2::div(uint8_t dividend, uint8_t divisor)
{
    if(divisor == 0) return 0; //illegal division by 0, but for now just return 0
    if(dividend == 0) return 0; //trivial division of 0
    //similarily to multiplication, x/y=b^(log(x)-log(y)), where b is the logarithm base
    return exp[(log[dividend] + 255 - log[divisor]) % 255];
}

/**
 * @brief Fast power in Galois field
 * @param x Base
 * @param exponent Exponent
 * @return Result
 */
uint8_t GF2::pow(uint8_t x, uint8_t exponent)
{
    //since a*log(x)=log(x^a) and b^log(x)=x, b^(a*log(x))=b^(log(x^a))=x^a, where b is the logarithm base
    return exp[(exponent * log[x]) % 255];
}

/**
 * @brief Fast inverse in Galois field
 * @param x Number of which inverse is calculated
 * @return 1/x
 */
uint8_t GF2::inv(uint8_t x)
{
    return exp[255 - log[x]];
}

/**
 * @brief Slow (no lookup table) multiplication algorithm in Galois field
 * @param x Multiplicand
 * @param y Multiplier
 * @return Multiplication result
 * Uses Russian Peasant Multiplication algorithm
 */
uint8_t GF2::slowMul(uint8_t x, uint8_t y)
{
    uint8_t ret = 0;
    uint16_t x_ = (uint16_t)x;
    while(y)
    {
        if(y & 1) //if current multiplier is odd
            ret ^= x_; //add multiplicand to the result
        y >>= 1; //divide y by 2
        x_ <<= 1; //multiply x_ by 2
        if(x_ & 256) x_ ^= GF2_POLY; //if there is a carry bit, apply modular reduction
    }
    return ret;
}

/**
 * @brief Check if object is initialized
 * @return 0 if initialized
 */
uint8_t isInitialized(void)
{
	return 0; //it must be initialized if the constructor was called
}


GF2::GF2()
{
    exp = new uint8_t[512];
    log = new uint8_t[256];

    uint8_t x = 1;
    //fill logarithm and exponential f. lookup tables for all possible values
    for(uint16_t i = 0; i < 256; i++)
    {
        exp[i] = x;
        log[x] = i;
        x = slowMul(x, 2);
    }
    for(uint16_t i = 256; i < 512; i++)
    {
        exp[i] = exp[i - 255]; //this is not necessary, but it will make things easier
    }
}

GF2::~GF2()
{
	if(exp != nullptr)
		delete [] exp;
	if(log != nullptr)
		delete [] log;
}
