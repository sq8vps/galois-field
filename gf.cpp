#include "gf.h"

/**
 * @brief Addition in Galois field
 * @param x Term 1
 * @param y Term 2
 * @return Sum
 */
uint16_t GF::add(uint16_t x, uint16_t y)
{
    return (x + y) % len;
}

/**
 * @brief Subtraction in Galois field
 * @param x Minuend
 * @param y Subtrahend
 * @return Difference
 */
uint16_t GF::sub(uint16_t x, uint16_t y)
{
    if(x >= y)
       return (x - y) % len;
    return (int32_t)len + (((int32_t)x - (int32_t)y) % (int32_t)len);
}

/**
 * @brief Multiplication in Galois field
 * @param x Multiplicand
 * @param y Multiplier
 * @return Multiplication result
 */
uint16_t GF::mul(uint16_t x, uint16_t y)
{
    if((x == 0) || (y == 0)) //trivial multiplication by 0
        return 0;
    return exp[log[x] + log[y]];
}

/**
 * @brief Division in Galois field
 * @param dividend Dividend
 * @param divisor Divisor
 * @return Division result. 0 is returned when dividing by 0.
 */
uint16_t GF::div(uint16_t dividend, uint16_t divisor)
{
    if(divisor == 0) return 0; //illegal division by 0, but for now just return 0
    if(dividend == 0) return 0; //trivial division of 0
    //similarly to multiplication, x/y=b^(log(x)-log(y)), where b is the logarithm base

    if((dividend % divisor) == 0)
        return ((dividend / divisor) % len);

    uint64_t a = dividend;
    for(uint16_t i = 1; i < divisor; i++)
    {
        a += len;
        if((a % divisor) == 0)
            return (a / divisor) % len;
    }
    return 0;
}

/**
 * @brief Power in Galois field
 * @param x Base
 * @param exponent Exponent
 * @return Result
 */
uint16_t GF::pow(uint16_t x, uint16_t exponent)
{
    //since a*log(x)=log(x^a) and b^log(x)=x, b^(a*log(x))=b^(log(x^a))=x^a, where b is the logairthm base
    return exp[(exponent * log[x]) % (len - 1)];
}

/**
 * @brief Inverse in Galois field
 * @param x Number of which inverse is calculated
 * @return 1/x
 */
uint16_t GF::inv(uint16_t x)
{
    return exp[(len - 1) - log[x]];
}

/**
 * @brief Slow (no lookup table) multiplication in Galois field
 * @param x Multiplicand
 * @param y Multiplier
 * @return Multiplication result
 */
uint16_t GF::slowMul(uint16_t x, uint16_t y)
{
    if(x == 0 || y == 0)
        return 0;

    return (x * y) % len;
}

#ifdef AAAAAAA

/**
 * @brief Multiplies polynomial by a scalar
 * @param p Input polynomial pointer
 * @param o Length of polynomial buffer (degree of a poly + 1)
 * @param s Scalar multiplicand
 * @param out Ouput polynomial. Has the same degree as input polynomial
 */
void ReedSolomon::gfPolyScale(uint8_t *p, uint16_t o, uint8_t s, uint8_t *out)
{
    for(uint16_t i = 0; i < o; i++)
    {
        out[i] = gfMul(p[i], s);
    }
}

/**
 * @brief Adds two polynomials
 * @param p1 1st polynomial
 * @param o1 1st polynomial buffer length (degree of a poly + 1)
 * @param p2 2nd polynomial
 * @param o2 2nd polynomial buffer length (degree of a poly + 1)
 * @param order Order of variables in polynomials (0 for c,x,x^2...)
 * @param out Output polynomial. Has the same degree as the input polynomial of a higher degree
 * @return Output polynomial length
 */
uint16_t ReedSolomon::gfPolyAdd(uint8_t *p1, uint16_t o1, uint8_t *p2, uint16_t o2, uint8_t *out, uint8_t order)
{
    if(o1 >= o2) //p1 is longer than p2
    {
        for(uint16_t i = 0; i < o1; i++)
        {
            out[i] = p1[i]; //copy longer (p1) polynomial
        }
        for(uint16_t i = (order ? (o1 - o2) : (0)); i < o2; i++)
        {
            out[i] ^= p2[i]; //sum (or subtract) the corresponding p2 coefficients
        }
    }
    else //p2 is longer than p1
    {
        for(uint16_t i = 0; i < o2; i++)
        {
            out[i] = p2[i];
        }
        for(uint16_t i = (order ? (o2 - o1) : (0)); i < o1; i++)
        {
            out[i] ^= p1[i];
        }
    }
    if(o1 >= o2) //the output polynomial has the same length as the longest input poly
        return o1;
    else
        return o2;
}

/**
 * @brief Multiplies two polynomials
 * @param p1 1st polynomial
 * @param o1 1st polynomial buffer length (degree of a poly + 1)
 * @param p2 2nd polynomial
 * @param o2 2nd polynomial buffer length (degree of a poly + 1)
 * @param out Output polynomial. Has a length of o1 + o2 - 1
 */
void ReedSolomon::gfPolyMul(uint8_t *p1, uint16_t o1, uint8_t *p2, uint16_t o2, uint8_t *out)
{
    memset(out, 0, o1 + o2 - 1);
    for(uint16_t i = 0; i < o1; i++)
    {
        for(uint16_t j = 0; j < o2; j++)
        {
            out[i+j] ^= gfMul(p1[i], p2[j]); //multiply each coefficient of p1 with each coefficient of p2 and sum the same powers
        }
    }
}

/**
 * @brief Evaluates the polynomial at the specified x
 * @param *p Polynomial
 * @param o Polynomial buffer length (degree of a poly + 1)
 * @param x Value to evaluate the polynomial at
 * @return Evaluated value
 */
uint8_t ReedSolomon::gfPolyEval(uint8_t *p, uint16_t o, uint8_t x)
{
    uint8_t ret = p[0];
    for(uint16_t i = 1; i < o; i++)
    {
        ret = gfAdd(gfMul(ret, x), p[i]); //this uses Horner's scheme to perform fast evaluation
    }
    return ret;
}


/**
 * @brief Divide two polynomials
 * @param p1 Divident polynomial
 * @param o1 Divident polynomial buffer length
 * @param p2 Divisor polynomial
 * @param o2 Divisor polynomial buffer length
 * @param q Quotient polynomial (length = divident length - divisor length + 1)
 * @param r Remainder polynomial (length = divisor length - 1)
 */
void ReedSolomon::gfPolyDiv(uint8_t *p1, uint16_t o1, uint8_t *p2, uint16_t o2, uint8_t *q, uint8_t *r)
{
    uint8_t *t = new uint8_t[o1];
    memcpy(t, p1, o1);
    for(uint16_t i = 0; i < (o1 - o2 + 1); i++)
    {
        if(p1[i] == 0) continue;

        for(uint16_t j = 1; j < o2; j++)
        {
            if(p2[j] == 0) continue;
            t[i + j] ^= gfMul(p2[j], p1[i]);
        }
    }
    memcpy(q, t, o1 - o2 + 1);
    memcpy(r, t + o1, o2 - 1);

    delete[] t;
}

/**
 * @brief Divide two polynomials (returing remainder only)
 * @param p1 Divident polynomial
 * @param o1 Divident polynomial buffer length
 * @param p2 Divisor polynomial
 * @param o2 Divisor polynomial buffer length
 * @param r Remainder polynomial (length = divisor length - 1)
 */
void ReedSolomon::gfPolyDiv(uint8_t *p1, uint16_t o1, uint8_t *p2, uint16_t o2, uint8_t *r)
{
    uint8_t *t = new uint8_t[o1];
    memcpy(t, p1, o1);
    for(uint16_t i = 0; i < (o1 - o2 + 1); i++)
    {
        uint8_t coef = t[i];
        if(coef == 0) continue;

        for(uint16_t j = 1; j < o2; j++)
        {
            if(p2[j] == 0) continue;
            t[i + j] ^= gfMul(p2[j], coef);
        }
    }
    memcpy(r, t + o1 - o2 + 1, o2 - 1);

    delete[] t;
}

/**
 * @brief Generates generator (encoding) polynomial
 * @param t Number of redundancy symbols (n-k)
 * @param out Output polynomial (length = t + 1)
 */
void ReedSolomon::gfGeneratorPoly(uint16_t t, uint8_t *out)
{
    if(t == 0)
        return;
    memset(out, 0, t + 1);
    out[0] = 1;
    uint8_t *tmp = new uint8_t[t + 1];
    memset(tmp, 0, t + 1);
    for(uint16_t i = 0; i < t; i++)
    {
        memcpy(tmp, out, i + 1);
        uint8_t h[2] = {1, gfPow(2, (uint8_t)i)}; //generate irreducible polynomial
        gfPolyMul(tmp, i + 1, h, 2, out);
    }
    delete[] tmp;
}

/**
 * @brief Reverses the order of variables in a polynomial (in-place)
 * @param p Polynomial
 * @param o Polynomial length (order + 1)
 */
void ReedSolomon::gfPolyInv(uint8_t *p, uint16_t o)
{
    for(uint16_t i = 0; i < (o >> 1); i++)
    {
        uint8_t tmp = p[i];
        p[i] = p[o - i - 1]; //perform in-place swap
        p[o - i - 1] = tmp;
    }
}
#endif
GF::GF(uint16_t p)
{
    len = p;

    exp = new uint16_t[len << 1]; //initialize tables for Galois field exponential and logarithmic functions lookup
    log = new uint16_t[len];

    uint16_t x = 1;
    //fill logarithm and exponential f. lookup tables for all possible values
    for(uint16_t i = 0; i < len; i++)
    {
        exp[i] = x;
        log[x] = i;
        x = slowMul(x, 16);
    }
    for(uint16_t i = len; i < (len << 1); i++)
    {
        exp[i] = exp[i - len - 1]; //this is not neccessary, but it will speed up the multiplication and division process by using only the lookup tables
    }



}

GF::~GF()
{
    if(exp != nullptr)
        delete[] exp;
    if(log != nullptr)
        delete[] log;
}
