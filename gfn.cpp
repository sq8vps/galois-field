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
* @file gfn.cpp
* @brief Simple Galois field library for GF(p)
* @version 1.1
* @author Piotr Wilkon <sq8vps@gmail.com>
* @copyright Copyright 2021 Piotr Wilkon, licensed under GNU GPLv3
**/

#include "gfn.h"

/**
 * @brief Addition in Galois field
 * @param x Term 1
 * @param y Term 2
 * @return Sum
 */
uint16_t GFn::add(uint16_t x, uint16_t y)
{
    return (x + y) % len; //the most trivial operation here. Just add and then keep the result in GF boundaries
}

/**
 * @brief Subtraction in Galois field
 * @param x Minuend
 * @param y Subtrahend
 * @return Difference
 */
uint16_t GFn::sub(uint16_t x, uint16_t y)
{
    if(x >= y)
       return (x - y) % len;
    //len+(x-y)=len-(y-x) to avoid some signed/unsigned casts
    return len - ((y - x) % len);
}

/**
 * @brief Multiplication in Galois field
 * @param x Multiplicand
 * @param y Multiplier
 * @return Multiplication result
 */
uint16_t GFn::mul(uint16_t x, uint16_t y)
{
    if((x == 0) || (y == 0)) //trivial multiplication by 0
        return 0;

    return exp[(log[x] + log[y]) % (len - 1)];
}

/**
 * @brief Division in Galois field
 * @param dividend Dividend
 * @param divisor Divisor
 * @return Division result. 0 is returned when dividing by 0.
 */
uint16_t GFn::div(uint16_t dividend, uint16_t divisor)
{
    if(divisor == 0) return 0; //illegal division by 0, but for now just return 0
    if(dividend == 0) return 0; //trivial division of 0
    //similarly to multiplication, x/y=b^(log(x)-log(y)), where b is the logarithm base

    int32_t t = log[dividend] - log[divisor]; //temporarily store

    if(t >= 0) //logarithm difference is positive
    	return exp[t];

    return exp[(len - 1) + t]; //logarithm difference is negative and so the index, can't use negative indexes, so convert it to positive
}

/**
 * @brief Power in Galois field
 * @param x Base
 * @param exponent Exponent
 * @return Result
 */
uint16_t GFn::pow(uint16_t x, uint16_t exponent)
{
    //since a*log(x)=log(x^a) and b^log(x)=x, b^(a*log(x))=b^(log(x^a))=x^a, where b is the logarithm base
    return exp[(exponent * log[x]) % (len - 1)];
}

/**
 * @brief Inverse in Galois field
 * @param x Number of which inverse is calculated
 * @return 1/x
 */
uint16_t GFn::inv(uint16_t x)
{
    if(x == 0) //0 has no inverse
    	return 0; //but return 0


    return exp[(len - 1) - log[x]]; //x^(-1)=b^(-log(x)), but we don't have negative indexes, so just start from the last value in table (which is at len-1)
}

/**
 * @brief Slow (no lookup table) multiplication in Galois field
 * @param x Multiplicand
 * @param y Multiplier
 * @return Multiplication result
 */
uint16_t GFn::slowMul(uint16_t x, uint16_t y)
{
    if(x == 0 || y == 0)
        return 0;

    return (x * y) % len;
}

/**
 * @brief Check if number is prime
 * @param x Input number
 * @return 0 if prime, -1 if not
 */
int8_t GFn::checkPrime(uint16_t x)
{
	if(x < 2)
		return -1; //definitely not primes

	for(uint16_t i = 2; i < (x >> 1); i++)
	{
		if((x % i) == 0) //divisible by something - not a prime
			return -1;
	}
	return 0; //probably prime
}

/**
 * @brief Finds the highest prime number, but smaller than max
 * @param max The limit
 * @return Prime number, 0 if fail
 */
uint16_t GFn::findPrime(uint16_t max)
{
	if(max < 2)
		return 0;
	if(max == 2)
		return 2;

	max--;
	for(; max > 0; max--)
	{
		if(checkPrime(max) == 0)
			return max;
	}

	return 0;
}


/**
 * @brief Check if object is initialized
 * @return 0 if initialized
 */
uint8_t GFn::isInitialized(void)
{
	if(len) //there is some characteristic set, so the object is initialized
		return 0;

	return 1;
}

/**
 * @brief Initializes Galois Field object
 * @param p Field characteristic GF(p), must be prime
 */
GFn::GFn(uint16_t p)
{
    len = 0;

	if(checkPrime(p) != 0)
    	return; //not a prime number

    len = p; //store characteristic

	//TODO: I don't know yet why, but to generator number to generate lookup tables must be the highest prime number lower than the GF characteristic
	//otherwise we will get non-unique values
    //in "standard" Galois fields GF(p^n), where n>1, the elements of this field are polynomials with a degree of up to n-1
    //the generator polynomial has a degree of n and must be irreducible
	uint16_t gen = findPrime(p);

    //initialize lookup tables for fast calculations
    exp = new uint16_t[len]; //exponential function table for every possible exponent in this field
    //if we have GF(p), there are p numbers in this field: 0,...,p-1
    //we can have the exp(x), where x is any element of GF(p), so it creates a table of p elements
    //for log(x) we can have all elements of GF(p) except 0 (log(0) is not defined).
    //the x (defined below) will wrap around to 1 in the last iteration
    //this means we have one non-unique value of x in our tables.
    //For example, in GF(7) with generator number 5, exp(0)=exp(6)=1 and that's true (7^0=1 and 7^6 mod 7=1)
    //although this is a problem for the logarithm, as it will have two different values for the same argument (log(1)=0 and log(1)=6)
    //log(x) is a function, so it must have only one value associated with one value. Just drop the log(1)=6.
    log = new uint16_t[len]; //logarithmic function table

    uint16_t x = 1;
    //fill lookup tables
    for(uint16_t i = 0; i < (len - 1); i++) //skip the last element for log table
    {
    	exp[i] = x;
        log[x] = i;
        x = slowMul(x, gen); //get next x by multiplying it by the generator number
    }
    exp[len - 1] = x; //store last element in exp table
}

GFn::~GFn()
{
    if(exp != nullptr)
        delete[] exp;
    if(log != nullptr)
        delete[] log;
}
