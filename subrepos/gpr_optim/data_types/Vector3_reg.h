/* ******************************************************************************
 * * MIT License *
 *                                                                                *
 * Copyright (c) 2018 Maxim Masterov *
 *                                                                                *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * * of this software and associated documentation files (the "Software"), to
 * deal  * in the Software without restriction, including without limitation the
 * rights   * to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell      * copies of the Software, and to permit persons to whom the
 * Software is          * furnished to do so, subject to the following
 * conditions:                       *
 *                                                                                *
 * The above copyright notice and this permission notice shall be included in
 * all * copies or substantial portions of the Software. *
 *                                                                                *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE  * SOFTWARE. *
 * ******************************************************************************
 */

#ifndef VECTOR3REG_H_
#define VECTOR3REG_H_

#include <smmintrin.h>
#include <x86intrin.h>

#include <cmath>
#include <iostream>

#ifdef __GNUG__
#ifndef _MM_ALIGN32
#define _MM_ALIGN32 __attribute__((aligned(32)))
#endif
#ifndef _MM_ALIGN16
#define _MM_ALIGN16 __attribute__((aligned(16)))
#endif
#define MUSTINLINE  //__attribute__((always_inline))
#define _aligned_malloc(x, y) malloc(x)
#define _aligned_free(x) free(x)
#else
#define MUSTINLINE __forceinline
#endif

// Define USE_FLOAT_VECTOR to use a single precision floating point format
// for representation of coordinates:
// #define USE_FLOAT_VECTOR

#ifdef USE_FLOAT_VECTOR
#define _ALIGN_SIZE _MM_ALIGN16
#else
#define _ALIGN_SIZE _MM_ALIGN32
#endif

namespace gpr {

/*!
 * @brief Regular class for 3d vector representation with double precision
 * floating points Class contains main methods for operating with 3d vectors.
 */
class _ALIGN_SIZE vector3_reg {
public:
#ifdef USE_FLOAT_VECTOR
    typedef float elt_type;
#else
    typedef double elt_type;
#endif

    union {
        struct _ALIGN_SIZE {
            elt_type x, y, z;
        };
        elt_type data[3];
    };

private:
    uint8_t been_inserted;  //!< Used in comma initializer

public:
    /*!
     * @brief default constructor
     */
    MUSTINLINE vector3_reg() : x(0), y(0), z(0), been_inserted(0) { }

    /*!
     * @brief Assign constructor
     */
    MUSTINLINE vector3_reg(elt_type _x, elt_type _y, elt_type _z)
        : x(_x), y(_y), z(_z), been_inserted(0)
    {
    }

    /*!
     * @brief Copy constructor
     */
    MUSTINLINE vector3_reg(const vector3_reg &other)
        : x(other.x), y(other.y), z(other.z), been_inserted(0)
    {
    }

    /*!
     * @brief Addition operator
     * @param other Other vector
     * @return Result of \e this type
     */
    MUSTINLINE vector3_reg operator+(const vector3_reg &other) const
    {
        return {x + other.x, y + other.y, z + other.z};
    }

    /*!
     * @brief Subtraction operator
     * @param other Other vector
     * @return Result of \e this type
     */
    MUSTINLINE vector3_reg operator-(const vector3_reg &other) const
    {
        return {x - other.x, y - other.y, z - other.z};
    }

    /*!
     * @brief Multiplication operator
     * @param other Other vector
     * @return Result of \e this type
     */
    MUSTINLINE vector3_reg operator*(const vector3_reg &other) const
    {
        return {x * other.x, y * other.y, z * other.z};
    }

    /*!
     * @brief Division operator
     * @param other Other vector
     * @return Result of \e this type
     */
    MUSTINLINE vector3_reg operator/(const vector3_reg &other) const
    {
        return {x / other.x, y / other.y, z / other.z};
    }

    /*!
     * @brief Addition assignment operator
     * @param other Other vector
     * @return Result of \e this type
     */
    MUSTINLINE vector3_reg &operator+=(const vector3_reg &other)
    {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    /*!
     * @brief Subtraction assignment operator
     * @param other Other vector
     * @return Result of \e this type
     */
    MUSTINLINE vector3_reg &operator-=(const vector3_reg &other)
    {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    /*!
     * @brief Multiplication assignment operator
     * @param other Other vector
     * @return Result of \e this type
     */
    MUSTINLINE vector3_reg &operator*=(const vector3_reg &other)
    {
        x *= other.x;
        y *= other.y;
        z *= other.z;
        return *this;
    }

    /*!
     * @brief Division assignment operator
     * @param other Other vector
     * @return Result of \e this type
     */
    MUSTINLINE vector3_reg &operator/=(const vector3_reg &other)
    {
        x /= other.x;
        y /= other.y;
        z /= other.z;
        return *this;
    }

    /*!
     * @brief Addition of scalar value to all coordinates of \e this vector
     * @param value value
     * @return Result of \e this type
     */
    MUSTINLINE vector3_reg operator+(elt_type value) const
    {
        return {x + value, y + value, z + value};
    }

    /*!
     * @brief Subtraction of scalar value to all coordinates of \e this vector
     * @param value value
     * @return Result of \e this type
     */
    MUSTINLINE vector3_reg operator-(elt_type value) const
    {
        return {x - value, y - value, z - value};
    }

    /*!
     * @brief Multiplication of all coordinates of \e this vector by scalar
     * value
     * @param value value
     * @return Result of \e this type
     */
    MUSTINLINE vector3_reg operator*(elt_type value) const
    {
        return {x * value, y * value, z * value};
    }
    friend MUSTINLINE vector3_reg operator*(elt_type value,
                                            const vector3_reg &rhs)
    {
        return rhs * value;
    }

    /*!
     * @brief Division of all coordinates of \e this vector by scalar value
     * @param value value
     * @return Result of \e this type
     */
    MUSTINLINE vector3_reg operator/(elt_type value) const
    {
        return {x / value, y / value, z / value};
    }

    /*!
     * @brief Addition and assignment operator for scalar value
     * Add given scalar value from all coordinates of \e this vector
     * @param value value
     * @return Reference to \e this vector
     */
    MUSTINLINE vector3_reg &operator+=(elt_type value)
    {
        x += value;
        y += value;
        z += value;
        return *this;
    }

    /*!
     * @brief Subtraction and assignment operator for scalar value
     * Subtract given scalar value from all coordinates of \e this vector
     * @param value value
     * @return Reference to \e this vector
     */
    MUSTINLINE vector3_reg &operator-=(elt_type value)
    {
        x -= value;
        y -= value;
        z -= value;
        return *this;
    }

    /*!
     * @brief Multiplication assignment operator for scalar value
     * Multiplies all coordinates of \e this vector by scalar value
     * @param value value
     * @return Reference to \e this vector
     */
    MUSTINLINE vector3_reg &operator*=(elt_type value)
    {
        x *= value;
        y *= value;
        z *= value;
        return *this;
    }

    /*!
     * @brief Division assignment operator for scalar value
     * Divides all coordinates of \e this vector by scalar value
     * @param value value
     * @return Reference to \e this vector
     */
    MUSTINLINE vector3_reg &operator/=(elt_type value)
    {
        x /= value;
        y /= value;
        z /= value;
        return *this;
    }

    /*!
     * @brief Assignment operator for scalars
     * Assigns given scalar to all coordinates
     * @param value Scalar value
     * @return Reference to \e this vector
     */
    MUSTINLINE vector3_reg &operator=(elt_type value)
    {
        x = value;
        y = value;
        z = value;
        return *this;
    }

    /*!
     * @brief Assignment operator for vector
     * Assigns given vector of three coordinates to \e this
     * @param other Other vector
     * @return Reference to \e this vector
     */
    MUSTINLINE vector3_reg &operator=(const vector3_reg &other)
    {
        x = other.x;
        y = other.y;
        z = other.z;
        return *this;
    }

    /*!
     * @brief Explicit set of three coordinates
     * Assigns given coordinates to \e this vector
     * @param _x Coordinate
     * @param _y Coordinate
     * @param _z Coordinate
     */
    MUSTINLINE void set(elt_type _x, elt_type _y, elt_type _z)
    {
        x = _x;
        y = _y;
        z = _z;
    }

    /*!
     * @brief Cross product of two vectors
     * @param other Other vector
     * @return Result of \e this type
     */
    MUSTINLINE vector3_reg cross(const vector3_reg &other) const
    {
        return {this->y * other.z - this->z * other.y,
                this->z * other.x - this->x * other.z,
                this->x * other.y - this->y * other.x};
    }

    /*!
     * @brief Dot product of two vectors
     * @param other Other vector
     * @return Result as a scalar
     */
    MUSTINLINE elt_type dot(const vector3_reg &other) const
    {
        return x * other.x + y * other.y + z * other.z;
    }

    /*!
     * @brief Length (absolute value) of \e this vector
     * @return Result as a scalar
     */
    MUSTINLINE elt_type length() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    /*!
     * @brief Reciprocal length (absolute value) of \e this vector
     * @return Result as a scalar
     */
    MUSTINLINE elt_type rlength() const
    {
        return 1. / length();
    }

    /*!
     * @brief Normalization of \e this vector
     * @return Vector scaled to unit length
     */
    MUSTINLINE vector3_reg normalize() const
    {
        vector3_reg v(x, y, z);
        elt_type abs = length();
        //      v = *this;
        v.x /= abs;
        v.y /= abs;
        v.z /= abs;
        return v;
    }

    /*!
     * @brief Sum of all elements
     * @return Sum of all elements
     */
    MUSTINLINE elt_type pack() const
    {
        return x + y + z;
    }

    /*!
     * @brief Prints coordinates of vector into the stream
     * @param os Reference to a stream
     * @param v Vector to be printed
     * @return Reference to the stream
     */
    friend MUSTINLINE std::ostream &operator<<(std::ostream &os,
                                               const vector3_reg &v)
    {
        os << v.x << ' ' << v.y << ' ' << v.z;  // << ' ';
        return os;
    }

    /*!
     * @brief Sets coordinates of vector from the stream
     * @param input Reference to a stream
     * @param v Vector to be used
     * @return Reference to the stream
     */
    friend MUSTINLINE std::istream &operator>>(std::istream &input,
                                               vector3_reg &v)
    {
        input >> v.x >> v.y >> v.z;
        return input;
    }

    /*!
     * @brief Allows to initialize vector in convenient way through operator<<
     */
    MUSTINLINE vector3_reg &operator<<(const elt_type &value)
    {
        x = value;
        been_inserted = 1;
        return *this;
    }

    /*!
     * @brief Inserts \e value in the vector
     * \warning There should be no more than 3 values passed into the vector
     * using comma initializer. Otherwise an error will be printed out
     */
    MUSTINLINE vector3_reg &operator,(const elt_type &value)
    {
        if (been_inserted == 0) {
            std::cerr << "Error! Bad use of comma initializer. Do not mess up "
                         "the code! "
                         "Assign values to vector using operator<<..."
                      << std::endl;
            return *this;
        }
        if (been_inserted == 1) {
            y = value;
            ++been_inserted;
        } else if (been_inserted == 2) {
            z = value;
            ++been_inserted;
        } else
            std::cerr << "Error! Too many arguments have been passed to comma "
                         "initializer "
                         "(see operator<<, vector3_reg)..."
                      << std::endl;
        return *this;
    }
};

} /* namespace gpr */

#endif /* VECTOR3REG_H_ */
