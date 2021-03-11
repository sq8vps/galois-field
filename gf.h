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
* @file gf.cpp
* @brief Simple Galois field library
* @author Piotr Wilkon <sq8vps@gmail.com>
* @copyright Copyright 2021 Piotr Wilkon, licensed under GNU GPLv3
**/


#ifndef GF_H
#define GF_H

#include <stdint.h>
#include <string.h>

class GF
{
public:
	/**
	 * @brief Addition in Galois field
	 * @param x Term 1
	 * @param y Term 2
	 * @return Sum
	 */
	uint16_t add(uint16_t x, uint16_t y);

	/**
	 * @brief Subtraction in Galois field
	 * @param x Minuend
	 * @param y Subtrahend
	 * @return Difference
	 */
	uint16_t sub(uint16_t x, uint16_t y);

	/**
	 * @brief Multiplication in Galois field
	 * @param x Multiplicand
	 * @param y Multiplier
	 * @return Multiplication result
	 */
	uint16_t mul(uint16_t x, uint16_t y);

	/**
	 * @brief Division in Galois field
	 * @param dividend Dividend
	 * @param divisor Divisor
	 * @return Division result. 0 is returned when dividing by 0.
	 */
	uint16_t div(uint16_t dividend, uint16_t divisor);

	/**
	 * @brief Power in Galois field
	 * @param x Base
	 * @param exponent Exponent
	 * @return Result
	 */
	uint16_t pow(uint16_t x, uint16_t exponent);

	/**
	 * @brief Inverse in Galois field
	 * @param x Number of which inverse is calculated
	 * @return 1/x
	 */
	uint16_t inv(uint16_t x);
	/**
	 * @brief Slow (no lookup table) multiplication in Galois field
	 * @param x Multiplicand
	 * @param y Multiplier
	 * @return Multiplication result
	 */
	uint16_t slowMul(uint16_t x, uint16_t y);

	/**
	 * @brief Check if number is prime
	 * @param x Input number
	 * @return 0 if prime, -1 if not
	 */
	static int8_t checkPrime(uint16_t x);

	/**
	 * @brief Finds the highest prime number, but smaller than max
	 * @param max The limit
	 * @return Prime number, 0 if fail
	 */
	static uint16_t findPrime(uint16_t max);
	/**
	 * @brief Initializes Galois Field object
	 * @param p Field characteristic GF(p), must be prime
	 */
	GF(uint16_t p);
	~GF();

private:
    uint16_t *exp; //exponent lookup table
    uint16_t *log; //logarithm lookup table
    uint16_t len; //field characteristic
};


#endif // REEDSOLOMON_H
