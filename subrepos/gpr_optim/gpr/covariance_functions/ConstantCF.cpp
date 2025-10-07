#include "ConstantCF.h"

namespace gpr {

ConstantCF::ConstantCF()
{
    constSigma2 = 0.;
}

double ConstantCF::calculateLogPrior()
{
    return 0;
}

Field<double> ConstantCF::calculateLogPriorGradient()
{
    //    Field<double>& lpg;
    return Field<double>();
}

void ConstantCF::calculateGradOfCovMatrix(const Coord &x, Coord &x2,
                                          std::vector<Field<double>> &DKff)
{
    Field<double> temp;

    if (x2.isEmpty()) x2 = x;
}

void ConstantCF::calculateGradOfCovMatrixWithDerivatives(
    const Coord &x, Coord &x2, Field<Index_t> &dims,
    std::vector<Field<double>> &DKff)
{
    Field<double> temp;
}

void ConstantCF::calculateGradOfCovMatrixWithDerivatives2(
    const Coord &x, Coord &x2, Field<Index_t> &dims1, Field<Index_t> &dims2,
    std::vector<Field<double>> &DKff)
{
    // This function does nothing!
    // See the MATLAB code, which has unsattisfied if-condition
}

void ConstantCF::ginput2(const Coord &x, Coord &x2, Field<Index_t> &dims,
                         std::vector<Field<double>> &DKff)
{
    Field<double> temp;

    DKff.resize(dims.getSize());
    temp.resize(x.getNumRows(), x2.getNumRows());

    for (Index_t i = 0; i < temp.getSize(); ++i)
        temp[i] = 0.;

    for (Index_t i = 0; i < dims.getSize(); ++i) {
        DKff[i].resize(temp.getNumRows(), temp.getNumCols());
        DKff[i] = temp;
    }
}

void ConstantCF::ginput3(const Coord &x, Coord &x2, Field<Index_t> &dims1,
                         Field<Index_t> &dims2,
                         std::vector<Field<double>> &DKff)
{
    Field<double> temp;

    DKff.resize(dims1.getSize() * dims2.getSize());
    temp.resize(x.getNumRows(), x2.getNumRows());

    for (Index_t i = 0; i < temp.getSize(); ++i)
        temp[i] = 0.;

    for (Index_t i = 0; i < (dims1.getSize() * dims2.getSize()); ++i) {
        DKff[i].resize(temp.getNumRows(), temp.getNumCols());
        DKff[i] = temp;
    }
}

void ConstantCF::ginput4(const Coord &x, Coord &x2, Field<Index_t> &dims,
                         std::vector<Field<double>> &DKff)
{
    Field<double> temp;

    DKff.resize(dims.getSize());
    temp.resize(x.getNumRows(), x2.getNumRows());

    for (Index_t i = 0; i < temp.getSize(); ++i)
        temp[i] = 0.;

    for (Index_t i = 0; i < dims.getSize(); ++i) {
        DKff[i].resize(temp.getNumRows(), temp.getNumCols());
        DKff[i] = temp;
    }
}

void ConstantCF::calculateCovarianceMatrix(const Coord &x, Coord &x2,
                                           Field<double> &C)
{
    Field<double> temp;
    temp.resize(x.getNumRows(), x2.getNumRows());
    C.resize(x.getNumRows(), x2.getNumRows());

    for (Index_t i = 0; i < temp.getSize(); ++i)
        temp[i] = this->constSigma2;

    C = temp;
}

} /* namespace gpr */
