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

#ifndef GF2_H
#define GF2_H

#include <stdint.h>

#define GF2_POLY 0x11d //primitive polynomial for division within the GF(2^8)

class GF2
{
public:
	/**
	 * @brief Addition in GF(2^8)
	 * @param x Term 1
	 * @param y Term 2
	 * @return Sum
	 */
	uint8_t add(uint8_t x, uint8_t y);

	/**
	 * @brief Subtraction in GF(2^8)
	 * @param x Minuend
	 * @param y Subtrahend
	 * @return Difference
	 */
	uint8_t sub(uint8_t x, uint8_t y);

	/**
	 * @brief Multiplication in GF(2^8)
	 * @param x Multiplicand
	 * @param y Multiplier
	 * @return Multiplication result
	 */
	uint8_t mul(uint8_t x, uint8_t y);

	/**
	 * @brief Division in GF(2^8)
	 * @param dividend Dividend
	 * @param divisor Divisor
	 * @return Division result. 0 is returned when dividing by 0.
	 */
	uint8_t div(uint8_t dividend, uint8_t divisor);

	/**
	 * @brief Power in GF(2^8)
	 * @param x Base
	 * @param exponent Exponent
	 * @return Result
	 */
	uint8_t pow(uint8_t x, uint8_t exponent);

	/**
	 * @brief Inverse in GF(2^8)
	 * @param x Number of which inverse is calculated
	 * @return 1/x
	 */
	uint8_t inv(uint8_t x);
	/**
	 * @brief Slow (no lookup table) multiplication in GF(2^8)
	 * @param x Multiplicand
	 * @param y Multiplier
	 * @return Multiplication result
	 */
	uint8_t slowMul(uint8_t x, uint8_t y);

	/**
	 * @brief Check if object is initialized
	 * @return 0 if initialized
	 */
	uint8_t isInitialized(void);

	/**
	 * @brief Initializes GF(2^8) object
	 */
	GF2();
	~GF2();

private:
    uint8_t *exp; //exponent lookup table
    uint8_t *log; //logarithm lookup table
};

#endif
