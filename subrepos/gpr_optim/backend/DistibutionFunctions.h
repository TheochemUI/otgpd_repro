//
//  DistibutionFunctions.h
//  gpr_dimer
//
//  Created by Maxim Masterov on 03/12/2020.
//

#ifndef DistibutionFunctions_h
#define DistibutionFunctions_h

#include <cmath>

namespace math {
/**
 * @brief Auxiliary functions for CDF calculation.
 */
class DistributionFunctions {
public:
    DistributionFunctions() { }
    ~DistributionFunctions() { }

    /**
     * @brief Normal cumulative distribution function.
     *
     * @param x Value for evaluation
     * @param mu Mean
     * @param sigma Standard deviation
     * @return Cumulative distribution function of the normal distribution
     */
    inline double normalCDF(const double x, double const mu, const double sigma)
    {
        double x_loc = (x - mu) / sigma;
        return 0.5 * erfc(-x_loc * M_SQRT1_2);
    }

    /**
     * @brief Inverse complementary error function
     *
     * @param p Probability value for evaluation
     * @param mu Mean
     * @param sigma Standard deviation
     * @return Inverse of the normal vumulative distribution function
     */
    inline double normalCDFInverse(const double p, double const mu,
                                   const double sigma)
    {
        return mu + sigma * M_SQRT2 * erfinv(2. * p - 1.);
    }

    /**
     * @brief Inverse error function
     * NOTE: this implementation is taken from
     * http://libit.sourceforge.net/math_8c-source.html
     * FIXME: check for other implementations
     * @param x Value for evaluation
     */
    double erfinv(double x);
};
}  // namespace math
#endif /* DistibutionFunctions_h */
