/*
 * Gradient.cpp
 *
 *  Created on: 20 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "Gradient.h"

namespace aux {

double Gradient::calculateDerivative(const double lengtScale2_rec,
                                     const double val1, const double val2)
{
    return -4. * lengtScale2_rec * val1 * val2;
}

void Gradient::resizeDerivativeFields(const gpr::Coord& x1,
                                      const gpr::Coord& x2,
                                      const gpr::AtomsConfiguration& conf_info,
                                      const uint8_t calc_options,
                                      const uint8_t field_type,
                                      gpr::Field<double>& D,
                                      std::vector<gpr::Field<double> >& D_pt)
{
    D.resize(x1.getNumRows(), x2.getNumRows());
    if (calc_options & field_type) {
        for (gpr::Index_t pt = 0; pt < conf_info.n_pt; ++pt)
            D_pt.push_back(D);
    }
}

void Gradient::evaluateValueInDerivativeField(
    const gpr::Pair<gpr::Index_t>& at, const uint8_t calc_options,
    const uint8_t field_type, const double lengtScale2_rec,
    const gpr::Index_t pt, const double val1, const double val2,
    gpr::Field<double>& D, std::vector<gpr::Field<double> >& D_pt)
{
    if (calc_options & field_type) {
        double deriv_ij = calculateDerivative(lengtScale2_rec, val1, val2);
        D(at.first, at.second) += deriv_ij;

        // Shift field type by 4 bits to the left to see if X_pt field
        // should be allocated
        if (calc_options & (field_type << 4))
            D_pt[pt](at.first, at.second) -= deriv_ij;
    }
}

void Gradient::pushBackDerivativeField(const uint8_t calc_options,
                                       const uint8_t field_type,
                                       gpr::Field<double>& D,
                                       std::vector<gpr::Field<double> >& D_vec)
{
    if (calc_options & field_type) D_vec.push_back(D);
}

void Gradient::calculateGradientBetweenMovingAtoms(
    const gpr::Coord& x1, const gpr::Coord& x2,
    const gpr::AtomsConfiguration& conf_info,
    const gpr::Field<double>& lengthScale, const gpr::Index_t direction_x1,
    const gpr::Index_t direction_x2, const gpr::Indices2D& ind,
    uint8_t calc_options, gpr::Derivatives<gpr::Field<double>*> derivatives,
    gpr::Derivatives<std::vector<gpr::Field<double> > >& derivatives_pt)
{
    gpr::Index_t n1 = x1.getNumRows();
    gpr::Index_t n2 = x2.getNumRows();

    gpr::Index_t pt = conf_info.pairtype(conf_info.atoms_mov.type(0, ind.i),
                                         conf_info.atoms_mov.type(0, ind.j));

    double lengtScale2_rec = 1. / (lengthScale[pt] * lengthScale[pt]);

    for (gpr::Index_t n = 0; n < n1; ++n) {
        double r_ij_1 = (x1.at(n, ind.j) - x1.at(n, ind.i)).length();
        double deriv_x1 = (x1.at(n, ind.i).data[direction_x1] -
                           x1.at(n, ind.j).data[direction_x1]) /
                          (r_ij_1 * r_ij_1 * r_ij_1);
        for (gpr::Index_t m = 0; m < n2; ++m) {
            double r_ij_2 = (x2.at(m, ind.j) - x2.at(m, ind.i)).length();
            double r_diff_rec = (1. / r_ij_1 - 1. / r_ij_2);
            double deriv_x2 = (x2.at(m, ind.i).data[direction_x2] -
                               x2.at(m, ind.j).data[direction_x2]) /
                              (r_ij_2 * r_ij_2 * r_ij_2);

            evaluateValueInDerivativeField(
                {n, m}, calc_options, OptionsForGradCalculation::D1,
                lengtScale2_rec, pt, deriv_x1, r_diff_rec, *derivatives.D1,
                derivatives_pt.D1);

            evaluateValueInDerivativeField(
                {n, m}, calc_options, OptionsForGradCalculation::D2,
                lengtScale2_rec, pt, -r_diff_rec, deriv_x2, *derivatives.D2,
                derivatives_pt.D2);

            evaluateValueInDerivativeField(
                {n, m}, calc_options, OptionsForGradCalculation::D12,
                lengtScale2_rec, pt, deriv_x1, deriv_x2, *derivatives.D12,
                derivatives_pt.D12);
        }
    }
}

void Gradient::calculateGradientBetweenMovingAndFrozenAtoms(
    const gpr::Coord& x1, const gpr::Coord& x2,
    const gpr::AtomsConfiguration& conf_info,
    const gpr::Field<double>& lengthScale, const gpr::Index_t direction_x1,
    const gpr::Index_t direction_x2, const gpr::Indices2D& ind,
    uint8_t calc_options, gpr::Derivatives<gpr::Field<double>*> derivatives,
    gpr::Derivatives<std::vector<gpr::Field<double> > >& derivatives_pt)
{
    gpr::Index_t n1 = x1.getNumRows();
    gpr::Index_t n2 = x2.getNumRows();

    gpr::Index_t pt =
        conf_info.pairtype(conf_info.atoms_mov.type(0, ind.i),
                           conf_info.atoms_froz_active.type(0, ind.j));

    double lengtScale2_rec = 1. / (lengthScale[pt] * lengthScale[pt]);

    gpr::vector3_reg froz_active_pos =
        conf_info.atoms_froz_active.positions.at(0, ind.j);
    double coord_x1 = froz_active_pos.data[direction_x1];
    double coord_x2 = froz_active_pos.data[direction_x2];
    for (gpr::Index_t n = 0; n < n1; ++n) {
        gpr::vector3_reg x1_vec = x1.at(n, ind.i);
        double r_ij_1 = (x1_vec - froz_active_pos).length();
        double r_ij_1_rec = 1. / r_ij_1;
        double deriv_x1 = (x1_vec.data[direction_x1] - coord_x1) * r_ij_1_rec *
                          r_ij_1_rec * r_ij_1_rec;

        for (gpr::Index_t m = 0; m < n2; ++m) {
            gpr::vector3_reg x2_vec = x2.at(m, ind.i);
            double r_ij_2 = (x2_vec - froz_active_pos).length();
            double r_ij_2_rec = 1. / r_ij_2;
            double r_diff_rec = (r_ij_1_rec - r_ij_2_rec);
            double deriv_x2 = (x2_vec.data[direction_x2] - coord_x2) *
                              r_ij_2_rec * r_ij_2_rec * r_ij_2_rec;

            evaluateValueInDerivativeField(
                {n, m}, calc_options, OptionsForGradCalculation::D1,
                lengtScale2_rec, pt, deriv_x1, r_diff_rec, *derivatives.D1,
                derivatives_pt.D1);

            evaluateValueInDerivativeField(
                {n, m}, calc_options, OptionsForGradCalculation::D2,
                lengtScale2_rec, pt, -r_diff_rec, deriv_x2, *derivatives.D2,
                derivatives_pt.D2);

            evaluateValueInDerivativeField(
                {n, m}, calc_options, OptionsForGradCalculation::D12,
                lengtScale2_rec, pt, deriv_x1, deriv_x2, *derivatives.D12,
                derivatives_pt.D12);
        }
    }
}

void Gradient::calculateDerivativesSameDim(
    const gpr::Coord& x1, const gpr::Coord& x2,
    const gpr::AtomsConfiguration& conf_info,
    const gpr::Field<double>& lengthScale, const gpr::Field<gpr::Index_t>& dims,
    uint8_t calc_options,
    gpr::Derivatives<std::vector<gpr::Field<double> > >& derivatives,
    gpr::Derivatives<std::vector<gpr::Field<double> > >& derivatives_pt)
{
    gpr::Field<double> D1;
    gpr::Field<double> D2;
    gpr::Field<double> D12;
    gpr::Index_t N_mov = conf_info.atoms_mov.type.getNumCols();
    gpr::Index_t N_fro = conf_info.atoms_froz_active.type.getNumCols();

    resizeDerivativeFields(x1, x2, conf_info, calc_options,
                           OptionsForGradCalculation::D1_pt, D1,
                           derivatives_pt.D1);

    resizeDerivativeFields(x1, x2, conf_info, calc_options,
                           OptionsForGradCalculation::D2_pt, D2,
                           derivatives_pt.D2);

    resizeDerivativeFields(x1, x2, conf_info, calc_options,
                           OptionsForGradCalculation::D12_pt, D12,
                           derivatives_pt.D12);

    for (gpr::Index_t dim_counter = 0; dim_counter < dims.getSize();
         ++dim_counter) {
        gpr::Index_t i_dim = dims[dim_counter];
        gpr::Index_t i = ceil(i_dim / 3.) - 1;
        gpr::Index_t xyz = i_dim - i * 3 - 1;

        // Gradient between moving atoms
        if (N_mov > 1) {
            for (gpr::Index_t j = 0; j < N_mov; ++j) {
                if (i != j) {
                    calculateGradientBetweenMovingAtoms(
                        x1, x2, conf_info, lengthScale, xyz, xyz, {i, j},
                        calc_options, {&D1, &D2, &D12}, derivatives_pt);
                }
            }
        }

        // Gradient from moving atoms to active frozen atoms
        if (N_fro > 0) {
            for (gpr::Index_t j = 0; j < N_fro; ++j) {
                calculateGradientBetweenMovingAndFrozenAtoms(
                    x1, x2, conf_info, lengthScale, xyz, xyz, {i, j},
                    calc_options, {&D1, &D2, &D12}, derivatives_pt);
            }
        }

        pushBackDerivativeField(calc_options, OptionsForGradCalculation::D1, D1,
                                derivatives.D1);

        pushBackDerivativeField(calc_options, OptionsForGradCalculation::D2, D2,
                                derivatives.D2);

        pushBackDerivativeField(calc_options, OptionsForGradCalculation::D12,
                                D12, derivatives.D12);
    }
}

void Gradient::calculateDerivativesDiffDim(
    const gpr::Coord& x1, const gpr::Coord& x2,
    const gpr::AtomsConfiguration& conf_info,
    const gpr::Field<double>& lengthScale,
    const gpr::Field<gpr::Index_t>& dims1,
    const gpr::Field<gpr::Index_t>& dims2, const uint8_t calc_options,
    const uint8_t direction_x2,
    gpr::Derivatives<std::vector<gpr::Field<double> > >& derivatives,
    gpr::Derivatives<std::vector<gpr::Field<double> > >& derivatives_pt)
{
    gpr::Field<double> D1;
    gpr::Field<double> D2;
    gpr::Field<double> D12;
    gpr::Field<double> s2;
    gpr::Index_t n1 = x1.getNumRows();
    gpr::Index_t n2 = x2.getNumRows();
    gpr::Index_t N_mov = conf_info.atoms_mov.type.getNumCols();
    gpr::Index_t N_fro = conf_info.atoms_froz_active.type.getNumCols();

    D1.resize(n1, n2);
    D2.resize(n1, n2);
    D12.resize(n1, n2);
    s2.resize(1, lengthScale.getSize());

    resizeDerivativeFields(x1, x2, conf_info, calc_options,
                           OptionsForGradCalculation::D1_pt, D1,
                           derivatives_pt.D1);

    resizeDerivativeFields(x1, x2, conf_info, calc_options,
                           OptionsForGradCalculation::D2_pt, D2,
                           derivatives_pt.D2);

    resizeDerivativeFields(x1, x2, conf_info, calc_options,
                           OptionsForGradCalculation::D12_pt, D12,
                           derivatives_pt.D12);

    for (gpr::Index_t n = 0; n < s2.getSize(); ++n)
        s2[n] = 1. / (lengthScale[n] * lengthScale[n]);

    for (gpr::Index_t dim_counter_1 = 0; dim_counter_1 < dims1.getSize();
         ++dim_counter_1) {
        for (gpr::Index_t dim_counter_2 = 0; dim_counter_2 < dims2.getSize();
             ++dim_counter_2) {
            gpr::Index_t i_dim_1 = dims1[dim_counter_1];
            gpr::Index_t i_1 = ceil(i_dim_1 / 3.) - 1;
            gpr::Index_t i_dim_2 = dims2[dim_counter_2];
            gpr::Index_t i_2 = ceil(i_dim_2 / 3.) - 1;
            gpr::Index_t xyz_1 = i_dim_1 - i_1 * 3 - 1;
            gpr::Index_t xyz_2 = i_dim_2 - i_2 * 3 - 1;
            gpr::Index_t i;

            if (i_1 != i_2) {
                i = i_1;
                if (N_mov > 1) {
                    for (gpr::Index_t j = 0; j < N_mov; ++j) {
                        if (i != j) {
                            calculateGradientBetweenMovingAtoms(
                                x1, x2, conf_info, lengthScale, xyz_1, xyz_1,
                                {i, j},
                                calc_options &
                                    (OptionsForGradCalculation::D1 |
                                     OptionsForGradCalculation::D1_pt),
                                {&D1, &D2, &D12}, derivatives_pt);
                        }
                    }
                }

                if (N_fro > 0) {
                    for (gpr::Index_t j = 0; j < N_fro; ++j) {
                        calculateGradientBetweenMovingAndFrozenAtoms(
                            x1, x2, conf_info, lengthScale, xyz_1, xyz_1,
                            {i, j},
                            calc_options & (OptionsForGradCalculation::D1 |
                                            OptionsForGradCalculation::D1_pt),
                            {&D1, &D2, &D12}, derivatives_pt);
                    }
                }

                i = i_2;
                if (N_mov > 1) {
                    for (gpr::Index_t j = 0; j < N_mov; ++j) {
                        if (i != j) {
                            calculateGradientBetweenMovingAtoms(
                                x1, x2, conf_info, lengthScale, xyz_2, xyz_2,
                                {i, j},
                                calc_options &
                                    (OptionsForGradCalculation::D2 |
                                     OptionsForGradCalculation::D2_pt),
                                {&D1, &D2, &D12}, derivatives_pt);
                        }
                    }
                }

                if (N_fro > 0) {
                    for (gpr::Index_t j = 0; j < N_fro; ++j) {
                        calculateGradientBetweenMovingAndFrozenAtoms(
                            x1, x2, conf_info, lengthScale, xyz_2, xyz_2,
                            {i, j},
                            calc_options & (OptionsForGradCalculation::D2 |
                                            OptionsForGradCalculation::D2_pt),
                            {&D1, &D2, &D12}, derivatives_pt);
                    }
                }

                calculateGradientBetweenMovingAtoms(
                    x1, x2, conf_info, lengthScale, xyz_1, xyz_2, {i_1, i_2},
                    calc_options & (OptionsForGradCalculation::D12 |
                                    OptionsForGradCalculation::D12_pt),
                    {&D1, &D2, &D12}, derivatives_pt);

                // We need to multiply D12 and D12_pt by -1, because
                // calculateGradientBetweenMovingAndFrozenAtoms function
                // evaluates gradient of x2 with opposite sign (according to the
                // Matlab code)
                for (gpr::Index_t n = 0; n < D12.getSize(); ++n) {
                    D12[n] *= -1;
                }

                if (calc_options & OptionsForGradCalculation::D12_pt) {
                    for (gpr::Index_t n = 0; n < derivatives_pt.D12.size();
                         ++n) {
                        for (gpr::Index_t m = 0;
                             m < derivatives_pt.D12[n].getSize(); ++m) {
                            derivatives_pt.D12[n][m] *= -1.;
                        }
                    }
                }
            } else {
                i = i_1;

                if (N_mov > 1) {
                    for (gpr::Index_t j = 0; j < N_mov; ++j) {
                        if (i != j) {
                            calculateGradientBetweenMovingAtoms(
                                x1, x2, conf_info, lengthScale, xyz_1, xyz_2,
                                {i, j}, calc_options, {&D1, &D2, &D12},
                                derivatives_pt);
                        }
                    }
                }

                if (N_fro > 0) {
                    for (gpr::Index_t j = 0; j < N_fro; ++j) {
                        calculateGradientBetweenMovingAndFrozenAtoms(
                            x1, x2, conf_info, lengthScale, xyz_1, xyz_2,
                            {i, j}, calc_options, {&D1, &D2, &D12},
                            derivatives_pt);
                    }
                }
            }

            pushBackDerivativeField(calc_options, OptionsForGradCalculation::D1,
                                    D1, derivatives.D1);

            pushBackDerivativeField(calc_options, OptionsForGradCalculation::D2,
                                    D2, derivatives.D2);

            pushBackDerivativeField(calc_options,
                                    OptionsForGradCalculation::D12, D12,
                                    derivatives.D12);
        }
    }
}

} /* namespace aux */
