/*
 * Gradient.h
 *
 *  Created on: 20 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_GRADIENT_H_
#define GPR_GRADIENT_H_

#include <vector>

#include "../../data_types/Coord.h"
#include "../../data_types/Field.h"
#include "../../structures/Structures.h"

namespace aux {

/**
 * @brief Methods to calculate gradients (derivatives) with respect to different
 * directions.
 */
class Gradient {
public:
    Gradient() { }
    virtual ~Gradient() { }

    /**
     * @brief Evaluate gradient of covariance function with respect to both
     * input variables x and x2 (in same dimension).
     *
     * The \e first member of \e derivatives will contain a vector of fields for
     * every derivative calculated for \e x (analog of D1). The \e second member
     * of \e derivatives will contain derivatives calculated for active
     * pairtypes.
     * @param derivatives Pair with the result.
     *
     * @param x1 Matrix of gpr::Coordinates
     * @param x2 Matrix of gpr::Coordinates
     * @param conf_info Configuration data for the GP model
     * @param lengthScale Length-scale
     * @param dims Dimensions for x1 and x2 (???)
     * @param calc_options Determines which derivative to calculate (D1, D2,
     * D12)
     * @param derivatives Field of derivatives calculated at each \e dims
     * @param derivatives_pt Field of pairtype derivatives
     */
    void calculateDerivativesSameDim(
        const gpr::Coord& x1, const gpr::Coord& x2,
        const gpr::AtomsConfiguration& conf_info,
        const gpr::Field<double>& lengthScale,
        const gpr::Field<gpr::Index_t>& dims, uint8_t calc_options,
        gpr::Derivatives<std::vector<gpr::Field<double> > >& derivatives,
        gpr::Derivatives<std::vector<gpr::Field<double> > >& derivatives_pt);

    /**
     * @brief Evaluate gradient of covariance function, of which has been taken
     * partial derivatives with respect to both input variables x and x2 with
     * respect to parameters (if specified).
     *
     * @param x1 Matrix of gpr::Coordinates
     * @param x2 Matrix of gpr::Coordinates
     * @param conf_info Configuration data for the GP model
     * @param lengthScale Length-scale
     * @param dims1 Dimensions for x1 (???)
     * @param dims2 Dimensions for x2 (???)
     * @param calc_options General type of the options to be calculated
     * @param direction_x2 Direction for the derivative calculation for \e x2
     * field.
     * @param derivatives Field of derivatives calculated at each \e dims
     * @param derivatives_pt Field of pairtype derivatives
     */
    void calculateDerivativesDiffDim(
        const gpr::Coord& x1, const gpr::Coord& x2,
        const gpr::AtomsConfiguration& conf_info,
        const gpr::Field<double>& lengthScale,
        const gpr::Field<gpr::Index_t>& dims1,
        const gpr::Field<gpr::Index_t>& dims2, const uint8_t calc_options,
        const uint8_t direction_x2,
        gpr::Derivatives<std::vector<gpr::Field<double> > >& derivatives,
        gpr::Derivatives<std::vector<gpr::Field<double> > >& derivatives_pt);

private:
    /**
     * @brief Calculate derivative
     *
     * This function calculates derivative in the following form:
     * \f[
     *      res = -4 / l^2 val1 val2
     * \f]
     * where \f l \f is the length-scale.
     *
     * @param lengtScale2_rec Value of the reciprocal squared length-scale
     * @param val1 First value
     * @param val2 Second value
     * @return
     */
    double calculateDerivative(const double lengtScale2_rec, const double val1,
                               const double val2);

    /**
     * @brief Evaluates the derivative between two values.
     *
     * This function adds the calculated derivative to the field \e D and
     * subtracts the derivative from the \e pt-th field of \e D_pt.
     *
     * @param at Row and column in the field of derivatives
     * @param calc_options General type of the options to be calculated
     * @param field_type Type of the field to be allocated (D1, D2, D12)
     * @param lengtScale2_rec Value of the reciprocal squared length-scale
     * @param pt Pairtype
     * @param val1 First value
     * @param val2 Second value
     * @param D Field of derivatives
     * @param D_pt Field of pairtype derivatives
     */
    void evaluateValueInDerivativeField(
        const gpr::Pair<gpr::Index_t>& at, const uint8_t calc_options,
        const uint8_t field_type, const double lengtScale2_rec,
        const gpr::Index_t pt, const double val1, const double val2,
        gpr::Field<double>& D, std::vector<gpr::Field<double> >& D_pt);

    /**
     * @brief Calculate gradient between moving atoms
     *
     * @param x1 Matrix of gpr::Coordinates
     * @param x2 Matrix of gpr::Coordinates
     * @param conf_info Configuration data for the GP model
     * @param lengthScale Length-scale
     * @param direction_x1 Direction of the derivative for the \e x1 field
     * @param direction_x2 Direction of the derivative for the \e x2 field
     * @param ind Matrix indices at which the result is stored
     * @param calc_options General type of the options to be calculated
     * @param derivatives Field of derivatives calculated at each \e dims
     * @param derivatives_pt Field of pairtype derivatives
     */
    void calculateGradientBetweenMovingAtoms(
        const gpr::Coord& x1, const gpr::Coord& x2,
        const gpr::AtomsConfiguration& conf_info,
        const gpr::Field<double>& lengthScale, const gpr::Index_t direction_x1,
        const gpr::Index_t direction_x2, const gpr::Indices2D& ind,
        uint8_t calc_options, gpr::Derivatives<gpr::Field<double>*> derivatives,
        gpr::Derivatives<std::vector<gpr::Field<double> > >& derivatives_pt);

    /**
     * @brief Calculate gradient between moving and frozen atoms
     *
     * @param x1 Matrix of gpr::Coordinates
     * @param x2 Matrix of gpr::Coordinates
     * @param conf_info Configuration data for the GP model
     * @param lengthScale Length-scale
     * @param direction_x1 Direction of the derivative for the \e x1 field
     * @param direction_x2 Direction of the derivative for the \e x2 field
     * @param ind Matrix indices at which the result is stored
     * @param calc_options General type of the options to be calculated
     * @param derivatives Field of derivatives calculated at each \e dims
     * @param derivatives_pt Field of pairtype derivatives
     */
    void calculateGradientBetweenMovingAndFrozenAtoms(
        const gpr::Coord& x1, const gpr::Coord& x2,
        const gpr::AtomsConfiguration& conf_info,
        const gpr::Field<double>& lengthScale, const gpr::Index_t direction_x1,
        const gpr::Index_t direction_x2, const gpr::Indices2D& ind,
        uint8_t calc_options, gpr::Derivatives<gpr::Field<double>*> derivatives,
        gpr::Derivatives<std::vector<gpr::Field<double> > >& derivatives_pt);

    /**
     * @brief Resize derivative field and filed of pairtype derivatives.
     *
     * @param x1 Matrix of gpr::Coordinates
     * @param x2 Matrix of gpr::Coordinates
     * @param conf_info Configuration data for the GP model
     * @param calc_options General type of the options to be calculated
     * @param field_type Type of the field to be allocated (D1, D2, D12)
     * @param D Allocated field of derivatives
     * @param D_pt Allocated field of pairtype derivatives
     */
    void resizeDerivativeFields(const gpr::Coord& x1, const gpr::Coord& x2,
                                const gpr::AtomsConfiguration& conf_info,
                                const uint8_t calc_options,
                                const uint8_t field_type, gpr::Field<double>& D,
                                std::vector<gpr::Field<double> >& D_pt);

    /**
     * @brief Append field \e D to vector of fields \e D_vec.
     *
     * @param calc_options General type of the options to be calculated
     * @param field_type Type of the field to be allocated (D1, D2, D12)
     * @param D Field of derivatives
     * @param D_vec Vector of fields of derivatives
     */
    void pushBackDerivativeField(const uint8_t calc_options,
                                 const uint8_t field_type,
                                 gpr::Field<double>& D,
                                 std::vector<gpr::Field<double> >& D_vec);
};

} /* namespace aux */

#endif /* GPR_GRADIENT_H_ */
