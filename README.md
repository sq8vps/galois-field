# Simple Galois field library in C++
A very simple Galois field library for GF(2^8) and GF(p)

## Short description

There are two separate libraries: the first one (gf2.h) is for handling typical GF(2^8) which is commonly used for standard Reed-Solomon coding.
The second one (gfn.h) if for GF(p), where p is any 16-bit prime number.

## Data representation

Galois (finite) field elements are represented as polynomial. For GF(2^8) there are polynomials of a degree lower than 8 and with coefficients being integers modulo 2 (they can have values either 0 or 1).
These are stored as 8-bit bytes, where every bit is a polynomial coefficient. The LSB represents a coefficient for x^0, while the MSB - for x^7. This is a very common approach.

GF(p) is in fact GF(p^1). This means there are polynomials of a degree 0 with coefficients being integers modulo p. Effectively, all elements are just numbers. 

## Capabilites

* Addition
* Subtraction
* Multiplication (using lookup tables or typical calculation)
* Division
* Power
* Inverse
* Prime number generation 

## Contributing

Any contributions are appreciated.

## Authors

* **Piotr Wilkoń** - *Initial work* - [sq8vps](https://github.com/sq8vps)

## License

This project is licensed under the GPLv3 License - see the LICENSE.md file for details
