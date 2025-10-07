/*
 * GaussianProcessRegression.inl
 *
 *  Created on: 17 Nov 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_ML_GAUSSIANPROCESSREGRESSION_INL_
#define GPR_ML_GAUSSIANPROCESSREGRESSION_INL_

namespace gpr {

inline Index_t GaussianProcessRegression::getNumberOfRepetitiveIndices(
    const Eigen::VectorXd& ind_Ddim)
{
    Index_t counter = 1;
    Index_t first_elt = ind_Ddim(0);

    for (Index_t n = 1; n < ind_Ddim.rows(); ++n) {
        if (ind_Ddim(n) == first_elt)
            ++counter;
        else
            break;
    }

    return counter;
}

template <typename CovFunc>
void GaussianProcessRegression::applyCovarianceFunction(
    const EigenMatrix& x, const Eigen::VectorXd& ind_Ddim,
    const Eigen::VectorXd& uDdim, CovFunc& cov_func, EigenMatrix& Ktemp)
{
    Field<double> Kff;
    std::vector<Field<double> > Kdf;
    std::vector<Field<double> > Kdf2;
    std::vector<Field<double> > D;
    Field<Index_t> dims1, dims2;
    Coord x_loc_0, x_loc_1, x_loc_2;
    Index_t number_of_rep_indices = getNumberOfRepetitiveIndices(ind_Ddim);

    x_loc_0.resize(number_of_rep_indices, (Index_t)x.cols() + 1);
    x_loc_1.resize(number_of_rep_indices, (Index_t)x.cols() + 1);
    x_loc_2.resize(number_of_rep_indices, (Index_t)x.cols() + 1);

    dims1.resize(1, 1);
    dims2.resize(1, 1);

    // One dimensional input
    if (x.cols() < 2) {
        // FIXME: this part is not not tested, because it is not executed in
        // CuH2 case!

        // FIXME: is it correct?
        dims1(0, 0) = 0;

        extractCoordinatesByIndex(x, ind_Ddim, 0, x_loc_0);
        extractCoordinatesByIndex(x, ind_Ddim, 1, x_loc_1);

        cov_func.calculateCovarianceMatrix(x_loc_0, x_loc_0, Kff);
        cov_func.ginput4(x_loc_1, x_loc_0, dims1, Kdf);
        cov_func.ginput2(x_loc_1, x_loc_1, dims1, D);

        // TODO: check if this allocation is correct
        // FIXME: should be done with Eigen
        // temp = [Kff Kdf{1}'; Kdf{1} D{1}];
        // assuming that:
        //      size(Kff,1) == size(Kdf{1},1)
        //      size(Kdf,1) == size(D{1},1)
        //      size(Kff Kdf{1}',2) == size(Kdf{1} D{1},2)
        //      size(Kdf,1) == size(Kdf,2)
        Ktemp.resize(Kff.getNumRows() + Kdf[0].getNumRows(),
                     Kff.getNumCols() + Kdf[0].getNumCols());
        Ktemp.setZero();

        for (Index_t i = 0; i < Ktemp.rows(); ++i) {
            for (Index_t j = 0; j < Ktemp.cols(); ++j) {
                // block by block
                if (i > 0 && i < Kff.getNumRows()) {
                    if (j > 0 && j < Kff.getNumCols()) {
                        Ktemp(i, j) = Kff(i, j);
                    }
                    if (j >= Kff.getNumCols() && j < Ktemp.cols()) {
                        Ktemp(i, j) = Kdf[0](j - Kff.getNumCols(), i);
                    }
                }

                if (i >= Kff.getNumRows() && i < Ktemp.rows()) {
                    if (j > 0 && j < Kdf[0].getNumCols()) {
                        Ktemp(i, j) = Kdf[0](i, j);
                    }
                    if (j >= Kdf[0].getNumCols() && j < Ktemp.cols()) {
                        Ktemp(i, j) = D[0](i, j - Kdf[0].getNumCols());
                    }
                }
            }
        }
    } else {
        Field<double> tmp_cov;

        Ktemp.resize(x.rows(), x.rows());
        Ktemp.setZero();

        extractCoordinatesByIndex(x, ind_Ddim, 0, x_loc_0);
        // the block of covariance matrix
        cov_func.calculateCovarianceMatrix(x_loc_0, x_loc_0, tmp_cov);

        // FIXME: pointless copy. Data types in `trcov` should be
        // changed to avoid this copy.
        assignBlockToMatrix(ind_Ddim, ind_Ddim, 0, 0, tmp_cov, Ktemp);

        for (Index_t n = 0; n < uDdim.size(); ++n) {
            dims1(0, 0) = uDdim[n];
            extractCoordinatesByIndex(x, ind_Ddim, uDdim[n], x_loc_1);

            // the blocks on the left side, below Kff
            cov_func.ginput4(x_loc_1, x_loc_0, dims1, Kdf);
            cov_func.ginput2(x_loc_1, x_loc_1, dims1, D);

            // Ktemp(ind_Ddim==uDdim(u1),ind_Ddim==0) = Kdf{1};
            assignBlockToMatrix(ind_Ddim, ind_Ddim, uDdim[n], 0, Kdf[0], Ktemp);

            // Ktemp(ind_Ddim==0,ind_Ddim==uDdim(u1)) = Kdf{1}';
            assignBlockToMatrix(ind_Ddim, ind_Ddim, 0, uDdim[n], Kdf[0], Ktemp,
                                true);

            // Ktemp(ind_Ddim==uDdim(u1),ind_Ddim==uDdim(u1)) = D{1};
            assignBlockToMatrix(ind_Ddim, ind_Ddim, uDdim[n], uDdim[n], D[0],
                                Ktemp);

            for (Index_t m = n + 1; m < uDdim.size(); ++m) {
                dims2(0, 0) = uDdim[m];
                extractCoordinatesByIndex(x, ind_Ddim, uDdim[m], x_loc_2);
                cov_func.ginput3(x_loc_1, x_loc_2, dims1, dims2, Kdf2);

                // Ktemp(ind_Ddim==uDdim(u1),ind_Ddim==uDdim2(u2)) = Kdf2{1};
                assignBlockToMatrix(ind_Ddim, ind_Ddim, uDdim[n], uDdim[m],
                                    Kdf2[0], Ktemp);

                // Ktemp(ind_Ddim==uDdim2(u2),ind_Ddim==uDdim(u1)) = Kdf2{1}';
                assignBlockToMatrix(ind_Ddim, ind_Ddim, uDdim[m], uDdim[n],
                                    Kdf2[0], Ktemp, true);
            }
        }
    }
}

template <typename CovFunc>
void GaussianProcessRegression::evaluateCovarianceFunction(
    const EigenMatrix& x1, const EigenMatrix& x2,
    const Eigen::VectorXd& ind_Ddim, const Eigen::VectorXd& uDdim,
    const Eigen::VectorXd& ind_Ddim2, const Eigen::VectorXd& uDdim2,
    CovFunc& cov_func, EigenMatrix& Ktemp)
{
    /* FIXME: For both sexpat and const functions the same matrices x1 and x2
     * and vectors ind_Ddim and ind_Ddim2 are passed. Therefore, the coordinate
     * extraction operation is redundant and in a loop. Heavy optimization is
     * possible here. Firstly with respect to memory creation and destruction
     * and secondly optimizing for cache hits later.
     *
     * NOTE: This is a preliminary working implementation. Need to optimize.
     */

    Coord x1_loc_0, x2_loc_0, x1_loc_1, x2_loc_1;
    Field<double> C_temp;
    Field<Index_t> dims1, dims2;
    std::vector<Field<double> > Kdf(1), Kdf2(1);

    Eigen::VectorXd vec1 = ind_Ddim, vec2 = ind_Ddim2;

    Ktemp.resize(x1.rows(), x2.rows());
    Ktemp.setZero();
    dims1.resize(1, 1);
    dims2.resize(1, 1);

    /* Map uses the same memory allocation as the vectors. */
    Eigen::Map<Eigen::ArrayXd> ind_Ddim_copy(vec1.data(), vec1.size()),
        ind_Ddim2_copy(vec2.data(), vec2.size());

    x1_loc_0.resize((Index_t)(ind_Ddim_copy == 0).count(),
                    (Index_t)x1.cols() + 1);
    x2_loc_0.resize((Index_t)(ind_Ddim2_copy == 0).count(),
                    (Index_t)x2.cols() + 1);

    extractCoordinatesByIndex(x1, ind_Ddim, 0, x1_loc_0);
    extractCoordinatesByIndex(x2, ind_Ddim2, 0, x2_loc_0);

    // % Non-derivative observation non-derivative prediction
    if ((ind_Ddim_copy == 0).count() > 0 && (ind_Ddim2_copy == 0).count() > 0) {
        cov_func.calculateCovarianceMatrix(x1_loc_0, x2_loc_0, C_temp);
        assignBlockToMatrix(ind_Ddim, ind_Ddim2, 0, 0, C_temp, Ktemp);
    }

    // % Non-derivative observation, derivative prediction
    if ((ind_Ddim_copy == 0).count() > 0) {
        for (Index_t i = 0; i < uDdim2.size(); ++i) {
            dims2(0, 0) = uDdim2[i];
            x2_loc_1.resize((Index_t)(ind_Ddim2_copy == uDdim2[i]).count(),
                            (Index_t)x2.cols() + 1);
            extractCoordinatesByIndex(x2, ind_Ddim2, uDdim2[i], x2_loc_1);

            // Kdf = gpcf.fh.ginput4(gpcf, x2(ind_Ddim2==uDdim2(u2),:),
            // x1(ind_Ddim==0,:), uDdim2(u2));
            cov_func.ginput4(x2_loc_1, x1_loc_0, dims2, Kdf);
            // Ktemp(ind_Ddim==0,ind_Ddim2==uDdim2(u2)) = Kdf{1}';
            assignBlockToMatrix(ind_Ddim, ind_Ddim2, 0, uDdim2[i], Kdf[0],
                                Ktemp, true);
        }
    }

    // % Derivative observation non-derivative prediction
    if ((ind_Ddim2_copy == 0).count() > 0) {
        for (Index_t i = 0; i < uDdim.size(); ++i) {
            dims1(0, 0) = uDdim[i];
            x1_loc_1.resize((Index_t)(ind_Ddim_copy == uDdim[i]).count(),
                            (Index_t)x1.cols() + 1);
            extractCoordinatesByIndex(x1, ind_Ddim, uDdim[i], x1_loc_1);

            // Kdf = gpcf.fh.ginput4(gpcf, x1(ind_Ddim==uDdim(u1),:),
            // x2(ind_Ddim2==0,:), uDdim(u1));
            cov_func.ginput4(x1_loc_1, x2_loc_0, dims1, Kdf);

            // Ktemp(ind_Ddim==uDdim(u1),ind_Ddim2==0) = Kdf{1};
            assignBlockToMatrix(ind_Ddim, ind_Ddim2, uDdim[i], 0, Kdf[0],
                                Ktemp);
        }
    }

    // % Derivative observation, derivative prediction
    for (Index_t i = 0; i < uDdim.size(); ++i) {
        dims1(0, 0) = uDdim[i];
        x1_loc_1.resize((Index_t)(ind_Ddim_copy == uDdim[i]).count(),
                        (Index_t)x1.cols() + 1);
        extractCoordinatesByIndex(x1, ind_Ddim, uDdim[i], x1_loc_1);
        for (Index_t j = 0; j < uDdim2.size(); ++j) {
            dims2(0, 0) = uDdim2[j];
            x2_loc_1.resize((Index_t)(ind_Ddim2_copy == uDdim2[j]).count(),
                            (Index_t)x2.cols() + 1);
            extractCoordinatesByIndex(x2, ind_Ddim2, uDdim2[i], x2_loc_1);
            if (uDdim[i] == uDdim2[j]) {
                // Kdf2 = gpcf.fh.ginput2(gpcf, x1(ind_Ddim==uDdim(u1),:)
                // ,x2(ind_Ddim2==uDdim2(u2),:), uDdim(u1));
                cov_func.ginput2(x1_loc_1, x2_loc_1, dims1, Kdf2);
            } else {
                // Kdf2 = gpcf.fh.ginput3(gpcf, x1(ind_Ddim==uDdim(u1),:)
                // ,x2(ind_Ddim2==uDdim2(u2),:), uDdim(u1), uDdim2(u2));
                cov_func.ginput3(x1_loc_1, x2_loc_1, dims1, dims2, Kdf2);
            }
            // Ktemp(ind_Ddim==uDdim(u1),ind_Ddim2==uDdim2(u2)) = Kdf2{1};
            assignBlockToMatrix(ind_Ddim, ind_Ddim2, uDdim[i], uDdim2[j],
                                Kdf2[0], Ktemp);
        }
    }
}

template <typename CovFunc>
void GaussianProcessRegression::calculateGradientWithCovFunc(
    const EigenMatrix& x, const Eigen::VectorXd& ind_Ddim,
    const Eigen::VectorXd& uDdim, const Eigen::VectorXd& b,
    const EigenMatrix& invC, CovFunc& cov_func, Field<double>& gdata,
    Field<double>& gprior)
{
    std::vector<Field<double> > DKffa;
    Coord x_loc_0, x_loc_1, x_loc_2;
    std::vector<EigenMatrix> DKffc;
    Field<Index_t> dims1, dims2;
    Field<double> gprior_cf;
    Index_t number_of_rep_indices = getNumberOfRepetitiveIndices(ind_Ddim);

    // FIXME: check, currently we assume that ind_Ddim has pair of indices
    x_loc_0.resize(number_of_rep_indices, (Index_t)x.cols() + 1);
    x_loc_1.resize(number_of_rep_indices, (Index_t)x.cols() + 1);
    x_loc_2.resize(number_of_rep_indices, (Index_t)x.cols() + 1);

    dims1.resize(1, 1);
    dims2.resize(1, 1);

    gprior_cf = cov_func.calculateLogPriorGradient();
    gprior_cf *= -1.;

    extractCoordinatesByIndex(x, ind_Ddim, 0, x_loc_0);
    cov_func.calculateGradOfCovMatrix(x_loc_0, x_loc_0, DKffa);

    if (DKffa.size() != 0) {
        gdata.resize((Index_t)DKffa.size());
        DKffc.resize(DKffa.size());
        for (Index_t elt = 0; elt < DKffa.size(); ++elt) {
            DKffc[elt].resize(x.rows(), x.rows());
            assignBlockToMatrix(ind_Ddim, ind_Ddim, 0, 0, DKffa[elt],
                                DKffc[elt]);
        }

        for (Index_t n = 0; n < uDdim.size(); ++n) {
            std::vector<Field<double> > Kdf;
            std::vector<Field<double> > Kdf2;
            std::vector<Field<double> > D;

            dims1(0, 0) = uDdim[n];
            extractCoordinatesByIndex(x, ind_Ddim, uDdim[n], x_loc_1);

            // the blocks on the left side, below Kff
            cov_func.calculateGradOfCovMatrixWithDerivatives(x_loc_1, x_loc_0,
                                                             dims1, Kdf);
            cov_func.calculateGradOfCovMatrixWithDerivatives2(x_loc_1, x_loc_1,
                                                              dims1, dims1, D);

            for (Index_t elt = 0; elt < DKffa.size(); ++elt) {
                // DKffc{inp}(ind_Ddim==uDdim(u1),ind_Ddim==0) = Kdf{inp};
                assignBlockToMatrix(ind_Ddim, ind_Ddim, uDdim[n], 0, Kdf[elt],
                                    DKffc[elt]);

                // DKffc{inp}(ind_Ddim==0,ind_Ddim==uDdim(u1)) = Kdf{inp}';
                assignBlockToMatrix(ind_Ddim, ind_Ddim, 0, uDdim[n], Kdf[elt],
                                    DKffc[elt], true);

                // DKffc{inp}(ind_Ddim==uDdim(u1),ind_Ddim==uDdim(u1)) = D{inp};
                assignBlockToMatrix(ind_Ddim, ind_Ddim, uDdim[n], uDdim[n],
                                    D[elt], DKffc[elt]);
            }

            for (Index_t m = n + 1; m < uDdim.size(); ++m) {
                dims2(0, 0) = uDdim[m];
                extractCoordinatesByIndex(x, ind_Ddim, uDdim[m], x_loc_2);
                cov_func.calculateGradOfCovMatrixWithDerivatives2(
                    x_loc_1, x_loc_2, dims1, dims2, Kdf2);

                for (Index_t elt = 0; elt < DKffa.size(); ++elt) {
                    // DKffc{inp}(ind_Ddim==uDdim(u1),ind_Ddim==uDdim2(u2)) =
                    // Kdf2{inp};
                    assignBlockToMatrix(ind_Ddim, ind_Ddim, uDdim[n], uDdim[m],
                                        Kdf2[elt], DKffc[elt]);

                    // DKffc{inp}(ind_Ddim==uDdim2(u2),ind_Ddim==uDdim(u1)) =
                    // Kdf2{inp}';
                    assignBlockToMatrix(ind_Ddim, ind_Ddim, uDdim[m], uDdim[n],
                                        Kdf2[elt], DKffc[elt], true);
                }
            }
        }

        // Evaluate the gradient with respect to covariance function parameters
        for (Index_t elt = 0; elt < DKffa.size(); ++elt) {
            double Bdl = b.dot(DKffc[elt] * b);
            double Cdl = 0.;

            for (Index_t i = 0; i < invC.rows(); ++i) {
                for (Index_t j = 0; j < invC.cols(); ++j) {
                    Cdl += invC(i, j) * DKffc[elt](i, j);
                }
            }
            gdata[elt] = 0.5 * (Cdl - Bdl);
        }
        gprior.append(gprior_cf);
    }
}

inline void GaussianProcessRegression::setJitterSigma2(const Index_t value)
{
    jitter_sigma2 = value;
}

inline void GaussianProcessRegression::setParameters(const GPRSetup& parameters)
{
    sigma2 = parameters.sigma2;
    jitter_sigma2 = parameters.jitter_sigma2;
    optimization_alg = parameters.optimization_alg;
}

inline LikGaussian* GaussianProcessRegression::getLikGaussian()
{
    return lik_gaussian;
}

inline SexpatCF* GaussianProcessRegression::getSexpAtCovarianceFunction()
{
    return sexpat_cov_func;
}

inline ConstantCF* GaussianProcessRegression::getConstantCovarianceFunction()
{
    return const_cov_fun;
}

inline Index_t GaussianProcessRegression::getNumberOfPotentialCalls()
{
    return num_of_potential_calls;
}

} /* namespace gpr */

#endif /* GPR_ML_GAUSSIANPROCESSREGRESSION_INL_ */
