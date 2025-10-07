/*
 * Coord.h
 *
 *  Created on: 6 Jul 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef DATA_TYPES_COORD_H_
#define DATA_TYPES_COORD_H_

#include "Field.h"
#include "Vector3_reg.h"

namespace gpr {
/**
 * @brief Wrapper over the \e Field class. Allows operating on 3D vectors stored
 * in a 1D std::vector.
 */
class Coord : public Field<double> {
    typedef vector3_reg elt_type;

public:
    Coord() : Field() { }
    virtual ~Coord() { }

    /**
     * @brief Return reference to the vector.
     * @param i Vector index.
     * @return Reference to the 3D vector.
     */
    inline elt_type at(const Index_t i) const
    {
        return getVector(3 * i);
    }

    /**
     * @brief Return reference to the vector.
     * @param i Row index.
     * @param j Column index.
     * @return Reference to the 3D vector.
     */
    inline elt_type at(const Index_t i, const Index_t j) const
    {
        return getVector(3 * j + this->_cols * i);
    }

    /**
     * @brief Assign coordinates of a vector to the corresponding location in
     *        the field.
     * @param i Vector index.
     * @param vec 3D vector.
     */
    inline void set(const Index_t i, const vector3_reg vec)
    {
        Index_t id = i * 3;
        operator[](id) = vec.x;
        operator[](id + 1) = vec.y;
        operator[](id + 2) = vec.z;
    }

    /**
     * @brief Assign coordinates of a vector to the corresponding location in
     *        the field.
     * @param i Row index.
     * @param j Column index.
     * @param vec 3D vector.
     */
    inline void set(const Index_t i, const Index_t j, const vector3_reg vec)
    {
        Index_t id = j * 3;
        operator()(i, id) = vec.x;
        operator()(i, id + 1) = vec.y;
        operator()(i, id + 2) = vec.z;
    }

    /**
     * @brief Returns number of 3D vectors in \e this field.
     * \note This function returns the number of 3D points, not number of
     * elements!
     * @return Number of 3D vectors.
     */
    inline Index_t getNumPoints() const
    {
        return (Index_t)(this->getSize() / 3);
    }

    /**
     * @brief Append 3D \e vector to \e this field.
     */
    inline void appendVector3D(const vector3_reg& vector)
    {
        Index_t old_nj = this->getNumCols();
        Coord tmp;
        Index_t ni = this->getNumRows() > 0 ? this->getNumRows() : 1;
        Index_t nj = this->getNumCols() + 3;

        tmp = *this;

        this->resize(ni, nj);

        // Copy old field as it was (block)
        for (Index_t i = 0; i < tmp.getNumRows(); ++i) {
            for (Index_t j = 0; j < tmp.getNumCols(); ++j) {
                this->operator()(i, j) = tmp(i, j);
            }
        }

        // Append the new data
        for (Index_t i = 0; i < ni; ++i) {
            set(i, old_nj / 3, vector);
        }

        tmp.clear();
    }

    /**
     * @brief Delete 3D vector from the column \e j.
     * \note This operation is expensive
     * @param j Column index.
     */
    inline void deleteColumnAt(const Index_t j)
    {
        this->deleteColumn(3 * j + 2);
        this->deleteColumn(3 * j + 1);
        this->deleteColumn(3 * j);
    }

    /**
     * @brief Copy elements of the provided \e std::vector to elements of \e
     * this field.
     * @param other Reference vector.
     */
    inline void operator=(const std::vector<double>& other)
    {
        Field<double>::operator=(other);
    }

    /**
     * @brief Copy elements of the \e other field to elements of \e this field.
     * \note \e this field will be resized according to the \e other field.
     * @param other Reference field.
     * @return Reference to \e this field.
     */
    inline Coord& operator=(const Field<double>& other)
    {
        Field<double>::operator=(other);
        return *this;
    }

    /**
     * @brief Assign the provided \e Eigen::Vector to \e this field.
     * @note \e this field will be resized according to the \e other field.
     * @param other Reference field.
     * @return Reference to \e this field.
     */
    inline Coord& operator=(const Eigen::VectorXd& other)
    {
        Field<double>::operator=(other);
        return *this;
    }

    /**
     * @brief Assign elements of the \e other field to \e this field.
     * \note \e this object will be resized according to the \e other object.
     * @param other Reference coordinates.
     * @return Reference to \e this coordinates.
     */
    inline Coord& operator=(const Coord& other)
    {
        this->allocationCopy(other);
        for (Index_t n = 0; n < getSize(); ++n)
            this->operator[](n) = other[n];

        return *this;
    }

    /**
     * @brief Normalize row \e i.
     * This function calculates the norm or row \e i and divides all elements in
     * the row by the calculated norm.
     * @param i Row index.
     */
    // FIXME: add unit test
    inline void normalizeRow(const Index_t i)
    {
        double norm = 0.;
        for (Index_t j = 0; j < this->getNumPoints(); ++j) {
            norm += this->at(i, j).dot(this->at(i, j));
        }
        norm = sqrt(norm);
        double norm_rec = 1. / norm;
        for (Index_t n = 0; n < this->getSize(); ++n) {
            this->operator[](n) *= norm_rec;
        }
    }

private:
    /**
     * @brief Returns a 3D vector from the specified index \e id.
     */
    inline elt_type getVector(Index_t id) const
    {
        return {this->data[id], this->data[id + 1], this->data[id + 2]};
    }
};

} /* namespace gpr */

#endif /* DATA_TYPES_COORD_H_ */
