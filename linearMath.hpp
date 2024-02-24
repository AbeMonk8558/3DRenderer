#pragma once

#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <immintrin.h>
#include <raylib.h>


template <typename T>
class Vec2;
template <typename T>
class Vec3;
template <typename T>
class Matrix44;


using Vec2f = Vec2<float>;
using Vec3f = Vec3<float>;
using Matrix44f = Matrix44<float>;


template <typename T>
class Vec2
{
public:
    T x, y;

    Vec2() : x(T(0)), y(T(0)) {}

    Vec2(T _x, T _y) : x(_x), y(_y) {}

    Vec2(const Vec2<T>& other) : x(other.x), y(other.y) {}

    Vec2(const Vector2& raylibVec) : x(static_cast<T>(raylibVec.x)), y(static_cast<T>(raylibVec.y)) {}

    Vec2<T> operator + (const Vec2<T>& right) const
    {
        return Vec2<T>(x + right.x, y + right.y);
    }

    Vec2<T> operator - (const Vec2<T>& right) const
    {
        return Vec2<T>(x - right.x, y - right.y);
    }

    Vec2<T> operator * (const T& scalar) const
    {
        return Vec2<T>(x * scalar, y * scalar);
    }

    Vec2<T> operator / (const T& scalar) const
    {
        return Vec2<T>(x / scalar, y / scalar);
    }

    Vec2<T> operator - () const
    {
        return Vec2<T>(-x, -y);
    }

    T length() const
    {
        return sqrt(x * x + y * y);
    }

    T lengthSq() const
    {
        return x * x + y * y;
    }

    T dot(const Vec2<T>& right) const
    {
        return x * right.x + y * right.y;
    }

    Vec2<T> normalize() const
    {
        T invLen = T(1) / length();

        return Vec2<T>(x * invLen, y * invLen);
    }

    T cross(const Vec2<T>& right) const
    {
        // Effectively the same as taking the determinant of a 2x2 matrix
        return x * right.y - y * right.x;
    }

    Vec2<T> surfaceNorm() const
    {
        Vec2<T> n(-y, x);
        if (this->cross(n) < 0)
            return -n;
            
        return n;
    }
};


template <typename T>
class Vec3
{
public:
    T x, y, z;

    Vec3() : x(T(0)), y(T(0)), z(T(0)) {}

    Vec3(T _x, T _y, T _z) : x(_x), y(_y), z(_z) {}

    Vec3(const Vec3<T>& other) : x(other.x), y(other.y), z(other.z) {}

    Vec3<T> operator + (const Vec3<T>& right) const
    {
        return Vec3<T>(x + right.x, y + right.y, z + right.z);
    }

    Vec3<T> operator - (const Vec3<T>& right) const
    {
        return Vec3<T>(x - right.x, y - right.y, z - right.z);
    }

    Vec3<T> operator * (const T& scalar) const
    {
        return Vec3<T>(x * scalar, y * scalar, z * scalar);
    }

    Vec3<T> operator / (const T& scalar) const
    {
        return Vec3<T>(x / scalar, y / scalar, z / scalar);
    }

    T length() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    T lengthSq() const
    {
        return x * x + y * y + z * z;
    }

    T dot(const Vec3<T>& right) const
    {
        return x * right.x + y * right.y + z * right.z;
    }

    Vec3<T> normalize() const
    {
        T invLen = T(1) / length();

        return Vec3<T>(x * invLen, y * invLen, z * invLen);
    }

    Vec3<T> cross(const Vec3<T>& right) const
    {
        return Vec3<T>(
            y * right.z - z * right.y,
            z * right.x - x * right.z,
            x * right.y - y * right.x
        );
    }
};


template <typename T>
class Matrix44
{
public:
    Matrix44() {}

    // Input values in a scalar-first fashion (e.g. do all x-values, then all y's, etc), 
    // rather than vector/axis first.
    Matrix44(std::initializer_list<T> values)
    {
        if (values.size() != 16)
            throw std::invalid_argument("A 4 by 4 matrix must be initialized with 16 values.");

        const T* v = values.begin();
        for (int i = 0; i < 16; i++)
        {
            _data[i / 4][i % 4] = *(v++);
        }
    }

    static Matrix44<T> identity()
    {
        return Matrix44<T>
        {
            T(1), T(0), T(0), T(0), 
            T(0), T(1), T(0), T(0), 
            T(0), T(0), T(1), T(0), 
            T(0), T(0), T(0), T(1)
        };
    }

    Matrix44<T> encode(const float& radians, const Vec3<T>& axis, const Vec3<T>& translation)
    {
        // TODO: implement this method
        return Matrix44<T>::identity();
    }

    const T* operator [] (const int& idx) const
    {
        return _data[idx];
    }

    T* operator [] (const int& idx)
    {
        return _data[idx];
    }

    Vec3<T> operator * (const Vec3<T>& v) const
    {
        T nx = v.x * _data[0][0] + v.y * _data[0][1] + v.z * _data[0][2] + _data[0][3];
        T ny = v.x * _data[1][0] + v.y * _data[1][1] + v.z * _data[1][2] + _data[1][3];
        T nz = v.x * _data[2][0] + v.y * _data[2][1] + v.z * _data[2][2] + _data[2][3];
        T nw = v.x * _data[3][0] + v.y * _data[3][1] + v.z * _data[3][2] + _data[3][3];

        if (nw != 1 && nw != 0)
        {
            nx /= nw;
            ny /= nw;
            nz /= nw;
        }

        return Vec3<T>(nx, ny, nz);
    }

    Matrix44<T> operator * (const Matrix44<T>& right) const
    {
        Matrix44<T> m;

        for (int r = 0; r < 4; r++)
        {
            for (int c = 0; c < 4; c++)
            {
                m[r][c] = _data[r][0] * right[0][c] + _data[r][1] * right[1][c] +
                        _data[r][2] * right[2][c] + _data[r][3] * right[3][c];
            }
        }

        return m;
    }

    Matrix44<T> transpose() const
    {
        Matrix44<T> transpose = *this;

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                transpose[j][i] = _data[i][j];

        return transpose;
    }

    // Uses Gauss-Jordan elimination and row elementary operations to calculate the inverse matrix.
    Matrix44<T> inverse() const
    {
        Matrix44<T> inv = identity();
        Matrix44<T> m = *this;

        for (int c = 0; c < 4; c++)
        {
            // Step 1: Make sure pivot coefficient is not 0, if so swap rows
            if (m[c][c] == T(0))
            {
                int max = c; // Pivot coefficient has row index equal to column index

                // Find row with maximum absolute value coefficient in same column
                for (int r = 0; r < 4; r++)
                    if (fabsf(m[r][c]) > fabsf(m[max][c])) max = r;

                if (max == c) return identity(); // TODO: Should probably throw exception or something

                for (int c2 = 0; c2 < 4; c2++)
                {
                    std::swap(m[max][c2], m[c][c2]);
                    std::swap(inv[max][c2], inv[c][c2]);
                }
            }

            // Step 2: Perform forward substitution on the column, setting all coefficients below
            //         the pivot within the column to 0 through row addition and scaling
            for (int r = c + 1; r < 4; r++) // Start at the row below the pivot in the column
            {
                T scalar = m[r][c] / m[c][c];

                for (int c2 = 0; c2 < 4; c2++)
                {
                    m[r][c2] -= scalar * m[c][c2];
                    inv[r][c2] -= scalar * inv[c][c2];
                }

                m[r][c] = 0; // Just to be safe
            }
        }

        // Step 3: Set all the pivot coefficients to 1 through row scaling
        for (int r = 0; r < 4; r++)
        {
            T scalar = T(1) / m[r][r];

            for (int c = 0; c < 4; c++)
            {
                m[r][c] *= scalar;
                inv[r][c] *= scalar;
            }

            m[r][r] = 1; // Just to be safe
        }

        // Step 4: Perform backward substitution on all columns above pivot
        for (int c = 0; c < 4; c++)
        {
            for (int r = 0; r < c; r++)
            {
                T scalar = m[r][c];

                for (int c2 = 0; c2 < 4; c2++)
                {
                    m[r][c2] -= scalar * m[c][c2];
                    inv[r][c2] -= scalar * inv[c][c2];
                }
            }
        }

        return inv;
    }

private:
    T _data[4][4] = 
    {
        {T(0), T(0), T(0), T(0)}, 
        {T(0), T(0), T(0), T(0)}, 
        {T(0), T(0), T(0), T(0)}, 
        {T(0), T(0), T(0), T(0)}
    };
};