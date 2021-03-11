#include <iostream>
#include "gf.h"

using namespace std;

int main()
{
    GF gf(1009);

    cout << gf.div(1, 3) << " " << gf.div(1, 1008);

    return 0;
}
