//
//  DistributionFunctions.cpp
//  gpr_dimer
//
//  Created by Maxim Masterov on 03/12/2020.
//

#include "DistibutionFunctions.h"

namespace math {

double DistributionFunctions::erfinv(double x)
{
    const double erfinv_a3 = -0.140543331;
    const double erfinv_a2 = 0.914624893;
    const double erfinv_a1 = -1.645349621;
    const double erfinv_a0 = 0.886226899;

    const double erfinv_b4 = 0.012229801;
    const double erfinv_b3 = -0.329097515;
    const double erfinv_b2 = 1.442710462;
    const double erfinv_b1 = -2.118377725;
    const double erfinv_b0 = 1;

    const double erfinv_c3 = 1.641345311;
    const double erfinv_c2 = 3.429567803;
    const double erfinv_c1 = -1.62490649;
    const double erfinv_c0 = -1.970840454;

    const double erfinv_d2 = 1.637067800;
    const double erfinv_d1 = 3.543889200;
    const double erfinv_d0 = 1;

    double x2, r, y;
    int sign_x;

    if (x < -1 || x > 1) return NAN;

    if (x == 0) return 0;

    if (x > 0)
        sign_x = 1;
    else {
        sign_x = -1;
        x = -x;
    }

    if (x <= 0.7) {
        x2 = x * x;
        r = x *
            (((erfinv_a3 * x2 + erfinv_a2) * x2 + erfinv_a1) * x2 + erfinv_a0);
        r /=
            (((erfinv_b4 * x2 + erfinv_b3) * x2 + erfinv_b2) * x2 + erfinv_b1) *
                x2 +
            erfinv_b0;
    } else {
        y = sqrt(-log((1 - x) / 2));
        r = (((erfinv_c3 * y + erfinv_c2) * y + erfinv_c1) * y + erfinv_c0);
        r /= ((erfinv_d2 * y + erfinv_d1) * y + erfinv_d0);
    }

    r = r * sign_x;
    x = x * sign_x;

    r -= (erf(r) - x) / (2 / sqrt(M_PI) * exp(-r * r));
    r -= (erf(r) - x) / (2 / sqrt(M_PI) * exp(-r * r));

    return r;
}
}  // namespace math
