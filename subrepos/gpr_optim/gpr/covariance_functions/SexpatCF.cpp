/*
 * SexpAt.cpp
 *
 *  Created on: 8 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "SexpatCF.h"

#include <cfloat>

#include "../auxiliary/Distance.h"
#include "../auxiliary/Gradient.h"
#include "../prior/PriorGaussian.h"
#include "../prior/PriorSqrtt.h"

namespace gpr {

SexpatCF::SexpatCF()
{
    clear();
}

SexpatCF::~SexpatCF()
{
    clear();
}

double SexpatCF::calculateLogPrior()
{
    PriorGaussian prior_gaussian;
    PriorSqrtt prior_sqrtt;
    double sum_length_scale = 0.;
    Field<double> magnSigma2_loc;
    double lp = 0.;

    prior_gaussian.setParameters(prior_parameters_gaussian);
    prior_sqrtt.setParameters(prior_parameters_sqrtt);

    // FIXME: check if magnSigma2 should be a Field, avoid unnecessary
    //        memory allocations
    magnSigma2_loc.resize(1, 1);
    magnSigma2_loc(0, 0) = magnSigma2;

    for (Index_t n = 0; n < lengthScale.getSize(); ++n) {
        sum_length_scale += log(lengthScale[n]);
    }

    lp += prior_sqrtt.calculateLogPrior(magnSigma2_loc) + log(magnSigma2);

    lp += prior_gaussian.calculateLogPrior(lengthScale) + sum_length_scale;

    return lp;
}

Field<double> SexpatCF::calculateLogPriorGradient()
{
    Field<double> lpg;
    PriorGaussian prior_gaussian;
    PriorSqrtt prior_sqrtt;
    Field<double> lpgs;
    Field<double> magnSigma2_loc;
    Index_t num_length_scales = lengthScale.getSize();

    prior_gaussian.setParameters(prior_parameters_gaussian);
    prior_sqrtt.setParameters(prior_parameters_sqrtt);

    // FIXME: check if magnSigma2 should be a Field, avoid unnecessary
    //        memory allocations
    lpg.resize(1, 1 + num_length_scales);
    magnSigma2_loc.resize(1, 1);
    magnSigma2_loc(0, 0) = magnSigma2;

    // FIXME: check if lpgs may have two elements
    lpgs = prior_sqrtt.calculateLogPriorGradient(magnSigma2_loc);
    lpg(0, 0) = lpgs[0] * magnSigma2 + 1.;
    //    lpg(0, 1) = lpgs[1];

    // FIXME: check if lpgs may have more then num_length_scales elements
    lpgs = prior_gaussian.calculateLogPriorGradient(lengthScale);
    for (Index_t n = 0; n < lengthScale.getSize(); ++n) {
        lpg(0, 1 + n) = lpgs[n] * lengthScale[n] + 1.;
        //    lpg(0, 2) = lpgs[num_length_scales+1];
    }

    return lpg;
}

void SexpatCF::calculateCovarianceMatrix(const Coord& x1, Coord& x2,
                                         Field<double>& C)
{
    // Based on equation (2.3)
    if (x2.isEmpty()) x2 = x1;

    aux::Distance dist_calculator;
    Field<double> dist_field;

    dist_calculator.dist_at(x1, x2, conf_info, lengthScale, dist_field);

    C.resize(dist_field.getNumRows(), dist_field.getNumCols());

    for (Index_t n = 0; n < C.getSize(); ++n) {
        // dist_field^2 because dist_at returns sqrt(dist)
        C[n] = magnSigma2 * exp(-dist_field[n] * dist_field[n] * 0.5);
        truncateValue(C[n]);
    }
}

void SexpatCF::calculateTrainingCovarianceMatrix(const Coord& x,
                                                 Field<double>& C)
{
    //    [n, m] =size(x);
    //    C = ones(n, 1).*gpcf.magnSigma2;
    //    C(C<eps)=0;

    // FIXME: bad allocation!
    C.resize(x.getNumRows(), 1);
    for (Index_t n = 0; n < C.getSize(); ++n) {
        C[n] = magnSigma2;
        truncateValue(C[n]);
    }
}

void SexpatCF::truncateValue(double& value)
{
    if (fabs(value) <= DBL_EPSILON) value = 0.;
}

void SexpatCF::calculateGradOfCovMatrix(const Coord& x1, Coord& x2,
                                        std::vector<Field<double> >& DKff)
{
    aux::Distance dist_calculator;
    io::ErrorManager err;
    Field<double> dist;
    std::vector<Field<double> > dist_pt;
    Field<double> K, DK_l;
    std::vector<Index_t> i1;

    DKff.clear();

    if (x2.isEmpty()) x2 = x1;

    // Evaluate: DKff{1} = d Kff / d magnSigma2
    //           DKff{2} = d Kff / d lengthScale
    // NOTE! Here we have already taken into account that the parameters
    // are transformed through log() and thus dK/dlog(p) = p * dK/dp
    // evaluate the gradient for training covariance
    if (x1.getNumCols() != x2.getNumCols()) {
        err << "gpcf_sexp -> _ghyper: The number of columns in x and x2 has to "
               "be the same.\n";
    }

    dist_calculator.dist_at(x1, x2, conf_info, lengthScale, dist);

    // We need squared distance here
    for (Index_t n = 0; n < dist.getSize(); ++n)
        dist[n] *= dist[n];

    K.resize(dist.getNumRows(), dist.getNumCols());

    for (Index_t n = 0; n < K.getSize(); ++n) {
        K[n] = magnSigma2 * exp(-0.5 * dist[n]);
        truncateValue(K[n]);
    }
    DKff.push_back(K);

    if (lengthScale.getSize() == 1) {
        // In the case of an isotropic sexp
        DK_l.resize(dist.getNumRows(), dist.getNumCols());
        for (Index_t n = 0; n < DK_l.getSize(); ++n)
            DK_l[n] = -magnSigma2 * exp(-0.5 * dist[n]) * dist[n];
        DKff.push_back(DK_l);
    } else {
        dist_calculator.dist_at_vec(x1, x2, conf_info, lengthScale, dist_pt);

        DK_l.resize(dist.getNumRows(), dist.getNumCols());
        for (Index_t m = 0; m < conf_info.n_pt; ++m) {
            for (Index_t n = 0; n < DK_l.getSize(); ++n) {
                DK_l[n] = -magnSigma2 * exp(-0.5 * dist[n]) * dist_pt[m][n];
            }
            DKff.push_back(DK_l);
        }
    }
}

void SexpatCF::calculateGradOfCovMatrixWithDerivatives(
    const Coord& x1, Coord& x2, Field<Index_t>& dims,
    std::vector<Field<double> >& DKff)
{
    std::vector<Field<double> > Cdm;
    aux::Distance dist_calculator;
    aux::Gradient grad_calculator;
    Field<double> dist;
    std::vector<Field<double> > dist_pt;
    std::vector<Field<double> > DK_pt;
    Derivatives<std::vector<Field<double> > > derivatives;
    Derivatives<std::vector<Field<double> > > derivatives_pt;
    Field<double> DK;
    uint8_t calc_options =
        OptionsForGradCalculation::D1 | OptionsForGradCalculation::D1_pt;

    DKff.clear();

    if (x2.isEmpty()) x2 = x1;

    ginput4(x1, x2, dims, Cdm);

    for (Index_t n = 0; n < Cdm.size(); ++n)
        DKff.push_back(Cdm[n]);

    dist_calculator.dist_at(x1, x2, conf_info, lengthScale, dist);
    dist_calculator.dist_at_vec(x1, x2, conf_info, lengthScale, dist_pt);

    grad_calculator.calculateDerivativesSameDim(x1, x2, conf_info, lengthScale,
                                                dims, calc_options, derivatives,
                                                derivatives_pt);

    // Note, `dist` and `D1` have the same size
    DK.resize(dist.getNumRows(), dist.getNumCols());

    for (Index_t n = 0; n < dims.getSize(); ++n) {
        for (Index_t pt = 0; pt < conf_info.n_pt; ++pt) {
            for (Index_t m = 0; m < DK.getSize(); ++m) {
                double expdist = exp(-0.5 * dist[m] * dist[m]);
                DK[m] = -magnSigma2 * expdist * derivatives_pt.D1[pt][m] +
                        0.5 * magnSigma2 * expdist * derivatives.D1[n][m] *
                            dist_pt[pt][m];
            }
            DK_pt.push_back(DK);
        }
    }

    for (Index_t pt = 0; pt < conf_info.n_pt; ++pt) {
        DKff.push_back(DK_pt[pt]);
    }
}

void SexpatCF::calculateGradOfCovMatrixWithDerivatives2(
    const Coord& x1, Coord& x2, Field<Index_t>& dims1, Field<Index_t>& dims2,
    std::vector<Field<double> >& DKff)
{
    std::vector<Field<double> > Cdm;
    aux::Distance dist_calculator;
    aux::Gradient grad_calculator;
    Field<double> dist;
    std::vector<Field<double> > dist_pt;
    std::vector<Field<double> > DK_pt;
    Derivatives<std::vector<Field<double> > > derivatives;
    Derivatives<std::vector<Field<double> > > derivatives_pt;
    std::vector<Field<double> > DKdd;
    Field<double> DK;
    uint8_t calc_options =
        OptionsForGradCalculation::D1 | OptionsForGradCalculation::D2 |
        OptionsForGradCalculation::D12 | OptionsForGradCalculation::D1_pt |
        OptionsForGradCalculation::D2_pt | OptionsForGradCalculation::D12_pt;

    DKff.clear();

    if (x2.isEmpty()) x2 = x1;

    if (dims1(0, 0) == dims2(0, 0))
        ginput2(x1, x2, dims1, DKdd);
    else
        ginput3(x1, x2, dims1, dims2, DKdd);

    DKff.push_back(DKdd[0]);

    dist_calculator.dist_at(x1, x2, conf_info, lengthScale, dist);
    dist_calculator.dist_at_vec(x1, x2, conf_info, lengthScale, dist_pt);

    grad_calculator.calculateDerivativesDiffDim(
        x1, x2, conf_info, lengthScale, dims1, dims2, calc_options,
        Directions::x, derivatives, derivatives_pt);

    DK.resize(dist.getNumRows(), dist.getNumCols());
    for (Index_t pt = 0; pt < conf_info.n_pt; ++pt) {
        for (Index_t n = 0; n < dims1.getSize(); ++n) {
            for (Index_t m = 0; m < dist.getSize(); ++m) {
                double expdist = exp(-0.5 * dist[m] * dist[m]);
                DK[m] = -magnSigma2 * expdist * derivatives_pt.D12[pt][m] -
                        0.25 * magnSigma2 * expdist * derivatives.D1[n][m] *
                            derivatives.D2[n][m] * dist_pt[pt][m];
                DK[m] += 0.5 * magnSigma2 * expdist *
                         (derivatives.D1[n][m] * derivatives_pt.D2[pt][m] +
                          derivatives.D2[n][m] * derivatives_pt.D1[pt][m] +
                          dist_pt[pt][m] * derivatives.D12[n][m]);
            }
        }
        DKff.push_back(DK);
    }
}

void SexpatCF::ginput2(const Coord& x1, Coord& x2, Field<Index_t>& dims,
                       std::vector<Field<double> >& DKff)
{
    aux::Distance dist_calculator;
    aux::Gradient grad_calculator;
    Field<double> dist;
    Derivatives<std::vector<Field<double> > > derivatives;
    Derivatives<std::vector<Field<double> > > derivatives_pt;
    Field<double> DK;
    Field<double> DK2;
    uint8_t calc_options = OptionsForGradCalculation::D1 |
                           OptionsForGradCalculation::D2 |
                           OptionsForGradCalculation::D12;

    DKff.clear();

    dist_calculator.dist_at(x1, x2, conf_info, lengthScale, dist);

    grad_calculator.calculateDerivativesSameDim(x1, x2, conf_info, lengthScale,
                                                dims, calc_options, derivatives,
                                                derivatives_pt);

    // Note, `dist` and `D1` have the same size
    DK.resize(dist.getNumRows(), dist.getNumCols());
    DK2.resize(dist.getNumRows(), dist.getNumCols());

    for (Index_t n = 0; n < dims.getSize(); ++n) {
        for (Index_t m = 0; m < dist.getSize(); ++m) {
            double expdist = exp(-0.5 * dist[m] * dist[m]);
            DK[m] = -0.5 * magnSigma2 * expdist * derivatives.D12[n][m];
            DK2[m] = 0.25 * magnSigma2 * expdist * derivatives.D1[n][m] *
                     derivatives.D2[n][m];
            DK[m] += DK2[m];
        }

        DKff.push_back(DK);
    }
}

void SexpatCF::ginput3(const Coord& x1, Coord& x2, Field<Index_t>& dims1,
                       Field<Index_t>& dims2, std::vector<Field<double> >& DKff)
{
    aux::Distance dist_calculator;
    aux::Gradient grad_calculator;
    Field<double> dist;
    Derivatives<std::vector<Field<double> > > derivatives;
    Derivatives<std::vector<Field<double> > > derivatives_pt;  // dummy object
    Field<double> DK;
    Field<double> DK2;
    uint8_t calc_options = OptionsForGradCalculation::D1 |
                           OptionsForGradCalculation::D2 |
                           OptionsForGradCalculation::D12;

    DKff.clear();

    dist_calculator.dist_at(x1, x2, conf_info, lengthScale, dist);

    grad_calculator.calculateDerivativesDiffDim(
        x1, x2, conf_info, lengthScale, dims1, dims2, calc_options,
        Directions::y, derivatives, derivatives_pt);

    // Note, `dist` and `D1` have the same size
    DK.resize(dist.getNumRows(), dist.getNumCols());
    DK2.resize(dist.getNumRows(), dist.getNumCols());

    Index_t n = 0;
    for (Index_t i = 0; i < dims1.getSize(); ++i) {
        for (Index_t j = 0; j < dims2.getSize(); ++j) {
            for (Index_t m = 0; m < dist.getSize(); ++m) {
                double expdist = exp(-0.5 * dist[m] * dist[m]);
                DK[m] = -0.5 * magnSigma2 * expdist * derivatives.D12[n][m];
                DK2[m] = 0.25 * magnSigma2 * expdist * derivatives.D1[n][m] *
                         derivatives.D2[n][m];
                DK[m] += DK2[m];
            }

            DKff.push_back(DK);
            ++n;
        }
    }
}

void SexpatCF::ginput4(const Coord& x1, Coord& x2, Field<Index_t>& dims,
                       std::vector<Field<double> >& DKff)
{
    aux::Distance dist_calculator;
    aux::Gradient grad_calculator;
    Field<double> dist;
    Derivatives<std::vector<Field<double> > > derivatives;
    Derivatives<std::vector<Field<double> > > derivatives_pt;
    Field<double> DK;
    uint8_t calc_options = OptionsForGradCalculation::D1;

    DKff.clear();

    if (x2.isEmpty()) x2 = x1;

    dist_calculator.dist_at(x1, x2, conf_info, lengthScale, dist);

    // We need squared distance here
    for (Index_t n = 0; n < dist.getSize(); ++n)
        dist[n] *= dist[n];

    grad_calculator.calculateDerivativesSameDim(x1, x2, conf_info, lengthScale,
                                                dims, calc_options, derivatives,
                                                derivatives_pt);

    // Note, `dist` and `D1` have the same size
    DK.resize(dist.getNumRows(), dist.getNumCols());

    for (Index_t n = 0; n < dims.getSize(); ++n) {
        for (Index_t m = 0; m < dist.getSize(); ++m) {
            DK[m] =
                -0.5 * magnSigma2 * exp(-0.5 * dist[m]) * derivatives.D1[n][m];
        }

        DKff.push_back(DK);
    }
}

} /* namespace gpr */
