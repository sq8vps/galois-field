# Simple Galois field library in C++
A very simple Galois field library for GF(p)

## Some desription

As said above, it handles Galois fields GF(p), where p is any 16-bit prime number. Doesn't support GF(p^n) (at least for now).
It was written to be used mainly in non-typical (not GF(2^8)) Reed-Solomon codes, but can be used in any solution.

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

This project is licensed under the GPLv3 License - see the [LICENSE.md](LICENSE.md) file for details
