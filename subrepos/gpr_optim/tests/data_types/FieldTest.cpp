/*
 * FieldTest.cpp
 *
 *  Created on: 26 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "FieldTest.h"

namespace gpr {
namespace tests {

FieldTest::FieldTest()
{
    // TODO Auto-generated constructor stub
}

FieldTest::~FieldTest()
{
    // TODO Auto-generated destructor stub
}

TEST_F(FieldTest, Resize1D)
{
    gpr::Field<double> field;
    gpr::Index_t size = 10;

    field.resize(size);

    EXPECT_EQ(field.getDimensions(), 1)
        << "The field has the wrong number of dimensions..";
    EXPECT_EQ(field.getSize(), size)
        << "The field has wrong number of elements.";
}

TEST_F(FieldTest, Resize2D)
{
    gpr::Field<double> field;
    gpr::Index_t size_i = 3;
    gpr::Index_t size_j = 10;

    field.resize(size_i, size_j);

    EXPECT_EQ(field.getDimensions(), 2)
        << "The field has the wrong number of dimensions..";
    EXPECT_EQ(field.getNumRows(), size_i)
        << "The field has wrong number of elements in i-th dimension.";
    EXPECT_EQ(field.getNumCols(), size_j)
        << "The field has wrong number of elements in j-st dimension.";
}

TEST_F(FieldTest, ResizeOneElement)
{
    gpr::Field<double> field;

    field.resize(1, 1, 1);

    EXPECT_EQ(field.getDimensions(), 1)
        << "The field has the wrong number of dimensions..";
    EXPECT_EQ(field.getSize(), 1) << "The field has wrong number of elements.";
}

TEST_F(FieldTest, ExtractDimension)
{
    gpr::Field<int> field1, field2;
    gpr::Index_t size_i = 3;
    gpr::Index_t size_j = 10;
    std::vector<int> ref(size_j);

    field1.resize(size_i, size_j);
    field2.resize(1, field1.getNumCols());

    ref = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    for (gpr::Index_t n = 0; n < field1.getNumRows(); ++n) {
        for (gpr::Index_t m = 0; m < field1.getNumCols(); ++m) {
            field1(n, m) = n + m;
        }
    }

    field2 = field1.extractRowAsVector(1);
    for (gpr::Index_t n = 0; n < field2.getNumCols(); ++n) {
        EXPECT_EQ(field2(0, n), ref[n]) << "The field has wrong values.";
    }
}

TEST_F(FieldTest, OperatorStarEqScalar)
{
    gpr::Field<int> field;
    gpr::Index_t size_i = 3;
    gpr::Index_t size_j = 10;
    std::vector<int> ref(size_j);
    int factor = -1;

    field.resize(size_i, size_j);

    ref = {-1, -2, -3, -4, -5, -6, -7, -8, -9, -10};

    for (gpr::Index_t n = 0; n < field.getNumRows(); ++n) {
        for (gpr::Index_t m = 0; m < field.getNumCols(); ++m) {
            field(n, m) = n + m;
        }
    }

    field *= factor;

    for (gpr::Index_t n = 0; n < field.getNumCols(); ++n) {
        EXPECT_EQ(field(1, n), ref[n]) << "The field has wrong values.";
    }
}

TEST_F(FieldTest, OperatorDivEqScalar)
{
    gpr::Field<double> field;
    gpr::Index_t size_i = 3;
    gpr::Index_t size_j = 10;
    std::vector<double> ref(size_j);
    double factor = -2;

    field.resize(size_i, size_j);

    ref = {-0.5, -1, -1.5, -2, -2.5, -3, -3.5, -4, -4.5, -5};

    for (gpr::Index_t n = 0; n < field.getNumRows(); ++n) {
        for (gpr::Index_t m = 0; m < field.getNumCols(); ++m) {
            field(n, m) = n + m;
        }
    }

    field /= factor;

    for (gpr::Index_t n = 0; n < field.getNumCols(); ++n) {
        EXPECT_EQ(field(1, n), ref[n]) << "The field has wrong values.";
    }
}

TEST_F(FieldTest, OperatorPlusEqVector)
{
    gpr::Field<int> field;
    gpr::Index_t size_i = 3;
    gpr::Index_t size_j = 10;
    std::vector<int> data(size_j);
    std::vector<int> ref(size_j);

    field.resize(size_i, size_j);

    data = {-2, -3, -4, -5, -6, -7, -8, -9, -10, -11};
    ref = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

    for (gpr::Index_t n = 0; n < field.getNumRows(); ++n) {
        for (gpr::Index_t m = 0; m < field.getNumCols(); ++m) {
            field(n, m) = n + m;
        }
    }

    field += data;

    for (gpr::Index_t n = 0; n < field.getNumCols(); ++n) {
        EXPECT_EQ(field(1, n), ref[n]) << "The field has wrong values.";
    }
}

TEST_F(FieldTest, OperatorMinusEqVector)
{
    gpr::Field<int> field;
    gpr::Index_t size_i = 3;
    gpr::Index_t size_j = 10;
    std::vector<int> data(size_j);
    std::vector<int> ref(size_j);

    field.resize(size_i, size_j);

    data = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    ref = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

    for (gpr::Index_t n = 0; n < field.getNumRows(); ++n) {
        for (gpr::Index_t m = 0; m < field.getNumCols(); ++m) {
            field(n, m) = n + m;
        }
    }

    field -= data;

    for (gpr::Index_t n = 0; n < field.getNumCols(); ++n) {
        EXPECT_EQ(field(1, n), ref[n]) << "The field has wrong values.";
    }
}

TEST_F(FieldTest, OperatorEqVector)
{
    gpr::Field<int> field;
    gpr::Index_t size_i = 3;
    gpr::Index_t size_j = 10;
    std::vector<int> data(size_j);

    field.resize(size_i, size_j);

    data = {-2, -3, -4, -5, -6, -7, -8, -9, -10, -11};

    field = data;

    for (gpr::Index_t n = 0; n < field.getNumCols(); ++n) {
        EXPECT_EQ(field(1, n), data[n]) << "The field has wrong values.";
    }
}

TEST_F(FieldTest, OperatorParentheses)
{
    gpr::Field<int> field;
    gpr::Index_t size_i = 3;
    gpr::Index_t size_j = 10;
    int ref_value = 5;

    field.resize(size_i, size_j);

    field(1, 6) = ref_value;

    EXPECT_EQ(field(1, 6), ref_value) << "The field has wrong value.";
}

TEST_F(FieldTest, AppendField)
{
    gpr::Field<int> field1, field2;
    gpr::Index_t size_i = 3;
    gpr::Index_t size_j = 4;
    std::vector<int> ref(2 * size_i * size_j);

    field1.resize(size_i, size_j);
    field2.resize(size_i, size_j);

    ref = {0, -1, -2, -3, 1, 0, -1, -2, 2, 1, 0, -1,
           0, 1,  2,  3,  1, 2, 3,  4,  2, 3, 4, 5};

    for (gpr::Index_t n = 0; n < field1.getNumRows(); ++n) {
        for (gpr::Index_t m = 0; m < field1.getNumCols(); ++m) {
            field1(n, m) = n + m;
            field2(n, m) = n - m;
        }
    }

    field2.append(field1);

    for (gpr::Index_t n = 0; n < field2.getSize(); ++n) {
        EXPECT_EQ(field2.getData()[n], ref[n]) << "The field has wrong values.";
    }
}

TEST_F(FieldTest, AppendFieldSlice)
{
    gpr::Field<int> field1, field2;
    gpr::Index_t size_i = 3;
    gpr::Index_t size_j = 4;
    std::vector<int> ref((1 + size_i) * size_j);

    field1.resize(size_i, size_j);
    field2.resize(size_i, size_j);

    ref = {0, -1, -2, -3, 1, 0, -1, -2, 2, 1, 0, -1, 2, 3, 4, 5};

    for (gpr::Index_t n = 0; n < field1.getNumRows(); ++n) {
        for (gpr::Index_t m = 0; m < field1.getNumCols(); ++m) {
            field1(n, m) = n + m;
            field2(n, m) = n - m;
        }
    }

    field2.append(field1, 2);

    for (gpr::Index_t n = 0; n < field2.getSize(); ++n) {
        EXPECT_EQ(field2.getData()[n], ref[n]) << "The field has wrong values.";
    }
}

TEST_F(FieldTest, Clear)
{
    gpr::Field<int> field1;
    gpr::Index_t size_i = 3;
    gpr::Index_t size_j = 4;

    field1.resize(size_i, size_j);

    for (gpr::Index_t n = 0; n < field1.getNumRows(); ++n) {
        for (gpr::Index_t m = 0; m < field1.getNumCols(); ++m) {
            field1(n, m) = n + m;
        }
    }

    field1.clear();

    EXPECT_EQ(field1.getSize(), 0)
        << "The field was not cleared (size is wrong).";
    EXPECT_EQ(field1.getNumRows(), 0)
        << "The field was not cleared (_ni is wrong).";
    EXPECT_EQ(field1.getNumCols(), 0)
        << "The field was not cleared (_nj is wrong).";
}

TEST_F(FieldTest, DeleteSlice)
{
    gpr::Field<int> field1;
    gpr::Index_t size_i = 3;
    gpr::Index_t size_j = 4;
    std::vector<int> ref((size_i - 1) * size_j);

    field1.resize(size_i, size_j);

    ref = {0, 1, 2, 3, 2, 3, 4, 5};

    for (gpr::Index_t n = 0; n < field1.getNumRows(); ++n) {
        for (gpr::Index_t m = 0; m < field1.getNumCols(); ++m) {
            field1(n, m) = n + m;
        }
    }

    field1.deleteRow(1);

    for (gpr::Index_t n = 0; n < field1.getSize(); ++n) {
        EXPECT_EQ(field1[n], ref[n]) << "The field has wrong values.";
    }
}

TEST_F(FieldTest, DeleteSliceV)
{
    gpr::Field<int> field1;
    gpr::Index_t size_i = 3;
    gpr::Index_t size_j = 4;
    std::vector<int> ref(size_i * (size_j - 1));

    field1.resize(size_i, size_j);

    ref = {0, 1, 3, 1, 2, 4, 2, 3, 5};

    for (gpr::Index_t n = 0; n < field1.getNumRows(); ++n) {
        for (gpr::Index_t m = 0; m < field1.getNumCols(); ++m) {
            field1(n, m) = n + m;
        }
    }

    field1.deleteColumn(2);

    for (gpr::Index_t n = 0; n < field1.getSize(); ++n) {
        EXPECT_EQ(field1[n], ref[n]) << "The field has wrong values.";
    }
}

} /* namespace tests */
} /* namespace gpr */
