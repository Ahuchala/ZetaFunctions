#include "flint/flint.h"
#include "flint/arb.h"

int main()
{
    // arb_t x;
    // arb_init(x);
    // arb_const_pi(x, 50 * 3.33);
    // arb_printn(x, 50, 0); flint_printf("\n");
    flint_printf("Computed with FLINT-%s\n", flint_version);
    // arb_clear(x);
}