/*
 * Field.h
 *
 *  Created on: 26 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef BACKEND_FIELD_H_
#define BACKEND_FIELD_H_

#include <algorithm>
#include <vector>

#include "../backend/Macros.h"
#include "../managers/io/ErrorManager.h"
#include "../managers/io/LogManager.h"
#include "Eigen/Dense"

namespace gpr {
/**
 * @brief Wrapper over the \e std::vector container. Allows operating on 2D
 * fields stored in a 1D vector. The data can be accessed using (i, j) indices.
 */
template <typename T>
class Field {
public:
    /**
     * @brief Deafult constructor.
     */
    Field() : _rows(0), _cols(0), _num_elts(0), _dimensions(0) { }

    /**
     * @brief Allocates memory for 1D \e Field and sets al values to 0.
     * @param rows Number of elements in i-th direction.
     */
    Field(const Index_t rows)
    {
        clear();
        resize(rows, 1, true);
    }

    /**
     * @brief Allocates memory for 2D \e Field and sets al values to 0.
     * @param rows Number of elements in i-th direction.
     * @param cols Number of elements in j-th direction.
     */
    Field(const Index_t rows, const Index_t cols)
    {
        clear();
        resize(rows, cols, true);
    }

    /**
     * @brief Copy constructor.
     */
    Field(const Field<T>& other)
    {
        clear();
        allocationCopy(other);

        for (Index_t n = 0; n < other.getSize(); ++n)
            data[n] = other[n];
    }

    /**
     * @brief Default distructor.
     */
    virtual ~Field()
    {
        clear();
    }

    /**
     * @brief Allocate memory for \e this field.
     * @param rows Number of elements in i-th dimension.
     * @param cols Number of elements in j-th dimension.
     * @param flush If true, set all elements to zero.
     */
    inline void resize(const Index_t rows, const Index_t cols = 1,
                       const bool flush = true)
    {
#ifndef NDEBUG
        assertMsg(
            rows != 0,
            "Number of elements in i-th direction cannot be equal to zero!");
        assertMsg(
            cols != 0,
            "Number of elements in j-th direction cannot be equal to zero!");
#endif

        _rows = rows;
        _cols = cols;
        calculateSize();
        _dimensions = countDimensions(rows, cols);

        data.resize(_rows * _cols);

        if (flush) setZero();
    }

    /**
     * @brief Return number of dimensions (1 or 2).
     */
    inline uint8_t getDimensions() const
    {
        return _dimensions;
    }

    /**
     * @brief Return total number of elements.
     */
    inline Index_t getSize() const
    {
        return (Index_t)data.size();
    }

    /**
     * @brief Return \e true if \e this field has no elements.
     */
    inline bool isEmpty() const
    {
        return getSize() == 0 ? true : false;
    }

    /**
     * @brief Return number of elements in i-th direction.
     */
    inline Index_t getNumRows() const
    {
        return _rows;
    }

    /**
     * @brief Return number of elements in j-th direction.
     */
    inline Index_t getNumCols() const
    {
        return _cols;
    }

    /**
     * @brief Return raw pointer to the data.
     */
    inline const T* getData() const
    {
        return data.data();
    }
    inline T* getData()
    {
        return data.data();
    }

    /**
     * @brief Return reference to the internal STL container.
     */
    inline const std::vector<T>& getInternalVector() const
    {
        return data;
    }
    inline std::vector<T>& getInternalVector()
    {
        return data;
    }

    /**
     * @brief Set all elements to zero.
     */
    inline void setZero()
    {
        for (Index_t n = 0; n < getSize(); ++n) {
            operator[](n) = 0;
        }
    }

    /**
     * @brief Set all elements in the specified row to zero.
     * @param row Index in the o-th direction.
     */
    inline void setZero(const Index_t row)
    {
        for (Index_t j = 0; j < getNumCols(); ++j) {
            operator()(row, j) = 0;
        }
    }

    /**
     * @brief Assign \e value to all elements of the field.
     * @param value Reference value.
     */
    inline void set(const T& value)
    {
        for (Index_t n = 0; n < getSize(); ++n) {
            operator[](n) = value;
        }
    }

    /**
     * @brief Return reference to the value using ID of the element.
     * @param id Element index.
     * @return Reference to the element.
     */
    inline const T& operator[](const Index_t id) const
    {
#ifndef NDEBUG
        assertMsg(id < getSize(), "Out of range! id = " + std::to_string(id) +
                                      ", size = " + std::to_string(getSize()));
#endif
        return data[id];
    }
    inline T& operator[](const Index_t id)
    {
        return const_cast<T&>(static_cast<const Field&>(*this).operator[](id));
    }

    /**
     * @brief Return reference to the value.
     * @param i Row index.
     * @param j Column index.
     * @return Reference to the element.
     */
    inline const T& operator()(const Index_t i, const Index_t j) const
    {
#ifndef NDEBUG
        assertMsg(i < _rows, "Out of range! i = " + std::to_string(i) +
                                 ", _ni = " + std::to_string(_rows));
        assertMsg(j < _cols, "Out of range! j = " + std::to_string(j) +
                                 ", _nj = " + std::to_string(_cols));
#endif
        return data[j + _cols * i];
    }
    inline T& operator()(const Index_t i, const Index_t j)
    {
        return const_cast<T&>(
            static_cast<const Field&>(*this).operator()(i, j));
    }

    /**
     * @brief Extract a row from \e this field.
     * @note This operation is inefficient!
     * @param i Row index.
     * @return Copy of elements from \e dimension.
     */
    inline const std::vector<T> extractRowAsVector(const Index_t i) const
    {
#ifndef NDEBUG
        assertMsg(i < _rows, "The specified data doesn't exist in the field!");
#endif

        std::vector<T> res(_cols);
        Index_t counter = _cols * i;

        for (Index_t j = 0; j < _cols; ++j) {
            res[j] = data[counter++];
        }

        return res;  // data[dimension];
    }

    /**
     * @brief Extract row as an \e Eigen::Vector.
     * @param i Row index.
     */
    inline Eigen::VectorXd extractRowAsEigenVector(const gpr::Index_t i) const
    {
        Eigen::VectorXd res;

        res.resize(getNumCols());

        for (Index_t n = 0; n < getNumCols(); ++n) {
            res(n) = this->operator()(i, n);
        }

        return res;
    }

    /**
     * @brief Assign elements of the provided vector to a row in \e this field.
     * @param i Row index.
     * @param other Reference vector.
     */
    inline void assignToRow(const Index_t i, const std::vector<T>& other)
    {
#ifndef NDEBUG
        assertMsg(i < _rows, "The specified data doesn't exist in the field!");
        assertMsg(
            _cols == other.size(),
            "The reference vector has more/less elements than this field!");
#endif

        Index_t counter = _cols * i;

        for (Index_t j = 0; j < _cols; ++j) {
            data[counter++] = other[j];
        }
    }

    /* ************************************************************ */
    /* ******************* Eigen based operations ***************** */
    /* ************************************************************ */

    /**
     * @brief Extract dense \e Eigen::Matrix matrix from \e this field.
     */
    inline Eigen::Matrix<T, EigenDynamic, EigenDynamic, EigenMatrixStorage>
    extractEigenMatrix() const
    {
        // Assume that Eigen supports RVO
        Eigen::Matrix<T, EigenDynamic, EigenDynamic, EigenMatrixStorage> result;

        result.resize(this->getNumRows(), this->getNumCols());

        for (Index_t i = 0; i < this->getNumRows(); ++i) {
            for (Index_t j = 0; j < this->getNumCols(); ++j) {
                result(i, j) = this->operator()(i, j);
            }
        }
        return result;
    }

    /**
     * @brief Extract dense \e Eigen::Vector vector from \e this field.
     */
    inline Eigen::Matrix<T, EigenDynamic, 1, EigenMatrixStorage>
    extractEigenVector() const
    {
        // Assume that Eigen supports RVO
        Eigen::Matrix<T, EigenDynamic, 1, EigenMatrixStorage> result;

        result.resize(this->getSize());

        for (Index_t n = 0; n < this->getSize(); ++n) {
            result(n) = this->operator[](n);
        }
        return result;
    }

    /**
     * @brief Assign the \e Eigen::Matrix to \this field.
     * @param reference_matrix Reference matrix.
     */
    inline void assignFromEigenMatrix(
        const Eigen::Matrix<T, EigenDynamic, EigenDynamic, EigenMatrixStorage>&
            reference_matrix)
    {
#ifndef NDEBUG
        assertMsg(this->getSize() == reference_matrix.size(),
                  "The size of Eigen matrix and Field do not match");
#endif
        for (Index_t i = 0; i < this->getNumRows(); ++i) {
            for (Index_t j = 0; j < this->getNumCols(); ++j) {
                this->operator()(i, j) = reference_matrix(i, j);
            }
        }
    }

    /* ************************************************************ */
    /* ******************* Assignment operators ******************* */
    /* ************************************************************ */
    /**
     * @brief Copy elements of the provided vector to elements of \e this field.
     * @param other Reference vector.
     */
    inline void operator=(const std::vector<T>& other)
    {
#ifndef NDEBUG
        assertMsg(
            _cols == other.size(),
            "The reference vector has more/less elements than this field!");
#endif

        for (Index_t i = 0; i < _rows; ++i) {
            for (Index_t j = 0; j < _cols; ++j) {
                operator()(i, j) = other[j];
            }
        }
    }

    /**
     * @brief Assign the provided field to \e this field.
     * @note \e this field will be resized according to the \e other field.
     * @param other Reference field.
     * @return Reference to \e this field.
     */
    inline Field<T>& operator=(const Field<T>& other)
    {
        if (other.getSize() == 0) {
            clear();
        } else {
            allocationCopy(other);
            for (Index_t n = 0; n < getSize(); ++n)
                operator[](n) = other[n];
        }
        return *this;
    }

    /**
     * @brief Assign the provided \e Eigen::Vector to \e this field.
     * @note \e this field will be resized according to the \e other field.
     * @param other Reference field.
     * @return Reference to \e this field.
     */
    inline Field<T>& operator=(const Eigen::VectorXd& other)
    {
        if (other.rows() == 0) {
            clear();
        } else {
            this->resize(1, (gpr::Index_t)other.rows());
            for (Index_t n = 0; n < getSize(); ++n)
                operator[](n) = other(n);
        }
        return *this;
    }

    /**
     * @brief Add elements of the provided \e std::vector to all elements of the
     * field. Addition is performed row-by-row.
     * @note Elements of the vector are added row-wise.
     * @note The size of the vector should be equivalent to the number of
     * elements in i-th direction in \e this field.
     * @return Reference to \e this field.
     */
    inline Field<T>& operator+=(const std::vector<T>& other)
    {
#ifndef NDEBUG
        assertMsg(_cols == other.size(),
                  "A mismatch between the size of the reference vector "
                  "and this field!");
#endif

        for (Index_t i = 0; i < _rows; ++i) {
            for (Index_t j = 0; j < _cols; ++j) {
                operator()(i, j) += other[j];
            }
        }
        return *this;
    }

    /**
     * @brief Subtract elements of the provided \e std::vector from all elements
     * of the field. Subtraction is performed row-by-row.
     * @note Elements of the vector are subtracted row-wise.
     * @note The size of the vector should be equivalent to the number of
     * elements in i-th direction in \e this field.
     * @return Reference to \e this field.
     */
    inline Field<T>& operator-=(const std::vector<T>& other)
    {
#ifndef NDEBUG
        assertMsg(_cols == other.size(),
                  "A mismatch between the size of the reference vector "
                  "and this field!");
#endif

        for (Index_t i = 0; i < _rows; ++i) {
            for (Index_t j = 0; j < _cols; ++j) {
                operator()(i, j) -= other[j];
            }
        }
        return *this;
    }

    /**
     * @brief Add elements of the \e other field to elements of \e this field.
     * @note Fields should have the same sizes.
     * @return Reference to \e this field.
     */
    inline Field<T>& operator+=(const Field<T>& other)
    {
#ifndef NDEBUG
        assertMsg(getSize() == other.getSize(),
                  "A mismatch between the size of the reference field "
                  "and this field!");
#endif

        for (Index_t n = 0; n < _num_elts; ++n) {
            operator[](n) += other[n];
        }
        return *this;
    }

    /**
     * @brief Subtract elements of the \e other field from elements of \e this
     * field.
     * @note Fields should have the same sizes.
     * @return Reference to \e this field.
     */
    inline Field<T>& operator-=(const Field<T>& other)
    {
#ifndef NDEBUG
        assertMsg(getSize() == other.getSize(),
                  "A mismatch between the size of the reference field "
                  "and this field!");
#endif

        for (Index_t n = 0; n < _num_elts; ++n) {
            operator[](n) -= other[n];
        }
        return *this;
    }

    /**
     * @brief Subtract elements of the \e other vector from elements of \e
     * this field.
     * @note Provided vector and \this field should have the same sizes.
     * @note The \e other vector can be of any type, as long as it implements
     * operator[] for element access.
     * @return Reference to \e this field.
     */
    template <typename D>
    inline Field<T>& operator-=(const D& other)
    {
        // #ifndef NDEBUG
        //         assertMsg(getSize() == other.getSize(),
        //                   "A mismatch between the size of the reference field
        //                   " "and this field!");
        // #endif

        for (Index_t n = 0; n < _num_elts; ++n) {
            operator[](n) -= other[n];
        }
        return *this;
    }

    /**
     * @brief Scale all elements of \e this field by \e value.
     * @param value Scaling factor.
     * @return Reference to \e this field.
     */
    inline Field<T>& operator*=(const T value)
    {
        for (Index_t n = 0; n < _num_elts; ++n) {
            operator[](n) *= value;
        }
        return *this;
    }

    /**
     * @brief Divide all elements of \e this field by \e value.
     * @param value Division factor.
     * @return Reference to \e this field.
     */
    inline Field<T>& operator/=(const T value)
    {
        return operator*=((T)1 / value);
    }

    /* ************************************************************ */
    /* ****************** Mathematical operators ****************** */
    /* ************************************************************ */
    /**
     * @brief Subtract elements of \e other field form elements of \e this field
     * and return the resullt.
     * @note Fields should have the same sizes.
     * @note \e this field remains intact.
     * @param other Reference field.
     */
    inline Field<T> operator-(const Field<T>& other) const
    {
#ifndef NDEBUG
        assertMsg(getSize() == other.getSize(),
                  "A mismatch between the size of the reference field "
                  "and this field!");
#endif

        Field<T> res;
        res.allocationCopy(*this);

        for (Index_t n = 0; n < _num_elts; ++n) {
            res[n] = operator[](n) - other[n];
        }

        return res;
    }

    /**
     * @brief Add elements of \e other field to elements of \e this field and
     * return the resullt.
     * @note Fields should have the same sizes.
     * @note \e this field remains intact.
     * @param other Reference field.
     */
    inline Field<T> operator+(const Field<T>& other) const
    {
#ifndef NDEBUG
        assertMsg(getSize() == other.getSize(),
                  "A mismatch between the size of the reference field "
                  "and this field!");
#endif

        Field<T> res;
        res.allocationCopy(*this);

        for (Index_t n = 0; n < _num_elts; ++n) {
            res[n] = operator[](n) + other[n];
        }

        return res;
    }

    /**
     * @brief Multiply elements of \e this field by elements of \e other field
     * and return the resullt.
     * @note Fields should have the same sizes.
     * @note \e this field remains intact.
     * @param other Reference field.
     */
    inline Field<T> operator*(const Field<T>& other) const
    {
#ifndef NDEBUG
        assertMsg(getSize() == other.getSize(),
                  "A mismatch between the size of the reference field "
                  "and this field!");
#endif

        Field<T> res;
        res.allocationCopy(*this);

        for (Index_t n = 0; n < _num_elts; ++n) {
            res[n] = operator[](n) * other[n];
        }

        return res;
    }

    /**
     * @brief Divide elements of \e this field by elements of \e other field and
     * return the resullt.
     * @note Fields should have the same sizes.
     * @note \e this field remains intact.
     * @param other Reference field.
     */
    inline Field<T> operator/(const Field<T>& other) const
    {
        // No need to duplicate this functionality
        return operator*((T)1 / other);
    }

    /**
     * @brief Subtract the given \e value from all elements of \e this field and
     * return the resullt.
     * @note \e this field remains intact.
     * @param value Reference value.
     * @return Copy of the result field.
     */
    inline Field<T> operator-(const T value) const
    {
        Field<T> res;

        res.allocationCopy(*this);

        for (Index_t n = 0; n < _num_elts; ++n) {
            res[n] = operator[](n) - value;
        }

        return res;
    }

    /**
     * @brief Add the given \e value from all elements of \e this field and
     * return the resullt and return the resullt.
     * @note \e this field remains intact.
     * @param value Reference value.
     * @return Copy of the result field.
     */
    inline Field<T> operator+(const T value) const
    {
        Field<T> res;

        res.allocationCopy(*this);

        for (Index_t n = 0; n < _num_elts; ++n) {
            res[n] = operator[](n) + value;
        }

        return res;
    }

    /**
     * @brief Scale \e this field and return the resullt.
     * @note \e this field remains intact.
     * @param value Scaling factor.
     * @return Copy of scaled field.
     */
    inline Field<T> operator*(const T value) const
    {
        Field<T> res;

        res.allocationCopy(*this);

        for (Index_t n = 0; n < _num_elts; ++n) {
            res[n] = operator[](n) * value;
        }

        return res;
    }

    /**
     * @brief Scale \e this field and return the resullt.
     * @note \e this field remains intact.
     * @param value Scaling factor.
     * @return Copy of scaled field.
     */
    inline Field<T> operator/(const T value) const
    {
        // No need to duplicate this functionality
        return operator*((T)1 / value);
    }

    /* ************************************************************ */
    /* ************** Other mathematical operations *************** */
    /* ************************************************************ */

    /**
     * @brief Calculate the dot product between two fields.
     * @note Fields should have the same sizes.
     * @param other Reference field.
     * @return Result of dot product.
     */
    inline T dot(const Field<T>& other) const
    {
#ifndef NDEBUG
        assertMsg(getSize() == other.getSize(),
                  "A mismatch between the size of the reference field "
                  "and this field!");
#endif

        T res;
        res = (T)0;
        for (Index_t n = 0, end = getSize(); n < end; ++n) {
            res += operator[](n) * other[n];
        }
        return res;
    }

    /**
     * @brief Calculate the dot product between two fields.
     * @note Fields should have the same sizes.
     * @param other Reference field.
     * @return Result of dot product.
     */
    inline T norm() const
    {
        T res;
        res = (T)0;
        for (Index_t n = 0, end = getSize(); n < end; ++n) {
            res += operator[](n) * operator[](n);
        }
        return std::sqrt(res);
    }

    /**
     * @brief Append \e other field to \e this field.
     * @note Fields should have the same dimensions and the same number of
     * columns.
     * @param other Reference field.
     */
    inline void append(const Field<T>& other)
    {
#ifndef NDEBUG
        assertMsg(other.getDimensions() != 0,
                  "Dimensions of the reference field cannot be zero!");
        if (!isEmpty())
            assertMsg(_cols == other.getNumCols(),
                      "The reference field has incompatible number of elements "
                      "in j-th "
                      "direction! _nj = " +
                          std::to_string(_cols) + ", other.getNumCols() = " +
                          std::to_string(other.getNumCols()));
#endif

        Index_t data_orig_size = (Index_t)data.size();

        if (isEmpty()) {
            _cols = other.getNumCols();
            countDimensions(_rows, _cols);
        }

        _rows += other.getNumRows();
        data.resize(data_orig_size + other.getSize());

        std::copy(other.getData(), other.getData() + other.getSize(),
                  data.begin() + data_orig_size);
    }

    /**
     * @brief Append a row \e i from the \e other field to \e this field.
     * @note Fields should have the same dimensions and the same number of
     * columns.
     * @param other Reference field.
     * @param i Row index of \e other field.
     */
    inline void append(const Field<T>& other, const Index_t i)
    {
#ifndef NDEBUG
        assertMsg(other.getDimensions() != 0,
                  "Dimensions of the reference field cannot be zero!");
        if (!isEmpty())
            assertMsg(_cols == other.getNumCols(),
                      "The reference field has incompatible number of elements "
                      "in j-th "
                      "direction! _nj = " +
                          std::to_string(_cols) + ", other.getNumCols() = " +
                          std::to_string(other.getNumCols()));
        assertMsg(i < other.getNumRows(),
                  "The specified index is greater than the field size! "
                  "i = " +
                      std::to_string(i) + ", other.getNumRows() = " +
                      std::to_string(other.getNumRows()));
#endif

        Index_t data_orig_size = (Index_t)data.size();
        Index_t start_pos_other = i * (Index_t)other.getNumCols();

        if (isEmpty()) {
            _cols = other.getNumCols();
            countDimensions(_rows, _cols);
        }

        incrementRows();

        data.resize(data_orig_size + other.getNumCols());

        std::copy(other.getData() + start_pos_other,
                  other.getData() + start_pos_other + other.getNumCols(),
                  data.begin() + data_orig_size);
    }

    /**
     * @brief Append \e other vector to \e this field.
     * @note The \e other vector can be \e std::vector or \e Eigen::VectorXd.
     * @param other Reference field.
     */
    template <typename D>
    inline void appendVector(const D& other)
    {
        Index_t data_orig_size = (Index_t)data.size();
        Index_t other_size = (Index_t)other.size();

        incrementRows();
        if (_cols == 0) _cols = other_size;
        if (_dimensions == 0) countDimensions(_rows, _cols);

        data.resize(data_orig_size + other_size);

        std::copy(other.data(), other.data() + other_size,
                  data.begin() + data_orig_size);
    }

    /**
     * @brief Append a row of \e value values to the end of \e this field in
     * i-th direction.
     * @param value Reference value.
     */
    inline void appendScalarAsRow(const T& value)
    {
        Index_t ni_orig = getNumRows();

        resize(ni_orig + 1, getNumCols(), 1, false);

        for (Index_t j = 0; j < getNumCols(); ++j)
            operator()(ni_orig, j) = value;
    }

    /**
     * @brief Delete indicated row from the field.
     * \warning This operation is computationally very expensive.
     * @param i Row index.
     */
    inline void deleteRow(const Index_t i)
    {
        if (isEmpty()) return;

#ifndef NDEBUG
        assertMsg(i < getNumRows(),
                  "The specified index is greater than the field size! "
                  "i = " +
                      std::to_string(i) +
                      ", getNumRows() = " + std::to_string(getNumRows()));
#endif
        data.erase(data.begin() + i * getNumCols(),
                   data.begin() + (i + 1) * getNumCols());
        decrementRows();

        // Flush all data if it was the last row
        if (getNumRows() == 0) clear();
    }

    /**
     * @brief Delete indicated column.
     * \warning This operation is computationally expensive.
     * @param j Column index.
     */
    inline void deleteColumn(const Index_t j)
    {
        if (isEmpty()) return;

#ifndef NDEBUG
        assertMsg(j < getNumCols(),
                  "The specified index is greater than the field size! "
                  "j = " +
                      std::to_string(j) +
                      ", getNumCols() = " + std::to_string(getNumCols()));
#endif
        Index_t pos = j + getNumCols() * (getNumRows() - 1);
        for (int n = pos; n >= 0; n -= getNumCols()) {
            data.erase(data.begin() + n);
        }
        decrementCols();

        // Flush all data if it was the last column
        if (getNumCols() == 0) clear();
    }

    /**
     * @brief Find maximum absolute value.
     */
    inline T getMaxAbsElt() const
    {
        T max_abs_val = (T)0;
        for (Index_t n = 0; n < getSize(); ++n) {
            T abs_val = std::fabs(operator[](n));
            if (abs_val > max_abs_val) {
                max_abs_val = abs_val;
            }
        }
        return max_abs_val;
    }

    /**
     * @brief Find maximum absolute value in a row.
     * @param i Row index.
     */
    inline T getMaxAbsEltInRow(const Index_t i) const
    {
#ifndef NDEBUG
        assertMsg(i < getNumRows(),
                  "The specified index is greater than the field size! "
                  "i = " +
                      std::to_string(i) +
                      ", getNumRows() = " + std::to_string(getNumRows()));
#endif

        T max_abs_val = (T)0;
        for (Index_t j = 0; j < getNumCols(); ++j) {
            T abs_val = std::fabs(operator()(i, j));
            if (abs_val > max_abs_val) {
                max_abs_val = abs_val;
            }
        }
        return max_abs_val;
    }

    /**
     * @brief Get mean value from the field.
     */
    inline T getMean() const
    {
        T res = (T)0;
        for (Index_t n = 0; n < getSize(); ++n) {
            res += operator[](n);
        }
        return res / getSize();
    }

    /**
     * @brief Calculate the Euclidean norm (magnitude) of the field.
     */
    inline T getNorm() const
    {
        T res = (T)0;
        for (Index_t n = 0; n < getSize(); ++n) {
            res += operator[](n) * operator[](n);
        }
        return std::sqrt(res);
    }

    /**
     * @brief Find maximum value in the field.
     */
    inline T getMaxElt() const
    {
        return *(std::max_element(data.begin(), data.end()));
    }

    /**
     * @brief Find minimum value in the field.
     */
    inline T getMinElt() const
    {
        return *(std::min_element(data.begin(), data.end()));
    }

    /**
     * @brief Clean up the field.
     */
    inline void clear()
    {
        data.clear();
        _rows = _cols = 0;
        _dimensions = 0;
        _num_elts = 0;
    }

    /**
     * @brief Print the field.
     */
    void print() const
    {
        io::LogManager log_man;

        log_man << "\n";
        for (Index_t i = 0; i < getNumRows(); ++i) {
            for (Index_t j = 0; j < getNumCols(); ++j) {
                log_man << operator()(i, j) << " ";
            }
            log_man << "\n";
        }
        log_man << "\n";
    }

    /**
     * @brief Print all sizes.
     */
    void printSizes() const
    {
        io::LogManager log_man;

        log_man << "getSize() = " << getSize() << "\n";
        log_man << "getNumRows() = " << getNumRows() << "\n";
        log_man << "getNumCols() = " << getNumCols() << "\n";
    }

protected:
    /**
     * @brief Copy allocation pattern from the provided field.
     * @param other Reference field.
     */
    inline void allocationCopy(const Field<T>& other)
    {
        if (other.isEmpty()) {
            clear();
            return;
        }
        resize(other.getNumRows(), other.getNumCols(), false);
    }

    /**
     * @brief Count number of dimensions.
     */
    inline uint8_t countDimensions(const Index_t ni, const Index_t nj)
    {
        uint8_t dimensions = 0;
        dimensions += ni != 1 ? 1 : 0;
        dimensions += nj != 1 ? 1 : 0;

        if (ni == nj && nj == 1) dimensions = 1;

        return dimensions;
    }

    /**
     * @ brief Recalculate Ni * Nj.
     */
    inline void calculateSize()
    {
        _num_elts = _rows * _cols;
    }

    /**
     * @brief Increment Ni.
     */
    inline void incrementRows()
    {
        ++_rows;
        calculateSize();
    }

    /**
     * @brief Decrement Ni.
     */
    inline void decrementRows()
    {
        --_rows;
        calculateSize();
    }

    /**
     * @brief Decrement Nj.
     */
    inline void decrementCols()
    {
        --_cols;
        calculateSize();
    }

protected:
    std::vector<T> data;  //!< Internal data.
    Index_t _rows;        //!< Number of elements in i-th direction.
    Index_t _cols;        //!< Number of elements in j-th direction.
    Index_t _num_elts;    //!< Multiplication of _ni and _nj.
    uint8_t _dimensions;  //!< Number of dimensions (1 or 2).
};

} /* namespace gpr */

#endif /* BACKEND_FIELD_H_ */
