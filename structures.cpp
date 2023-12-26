#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <raylib.h>
#include "3DRenderer.hpp"

// ************************** Vec2f ***************************

Vec2f::Vec2f() : x(0), y(0) {}

Vec2f::Vec2f(float _x, float _y) : x(_x), y(_y) {}

Vec2f::Vec2f(const Vec2f& other) : x(other.x), y(other.y) {}

Vec2f::Vec2f(const Vector2& raylibVec) : x(raylibVec.x), y(raylibVec.y) {}

Vec2f Vec2f::operator + (const Vec2f& right) const
{
    return Vec2f(x + right.x, y + right.y);
}

Vec2f Vec2f::operator - (const Vec2f& right) const
{
    return Vec2f(x - right.x, y - right.y);
}

Vec2f Vec2f::operator * (const float& scalar) const
{
    return Vec2f(x * scalar, y * scalar);
}

Vec2f Vec2f::operator / (const float& scalar) const
{
    return Vec2f(x / scalar, y / scalar);
}

float Vec2f::length() const
{
    return sqrtf(x * x + y * y);
}

float Vec2f::lengthSq() const
{
    return x * x + y * y;
}

float Vec2f::dot(const Vec2f& right) const
{
    return x * right.x + y * right.y;
}

Vec2f Vec2f::normalize() const
{
    float invLen = 1 / length();

    return Vec2f(x * invLen, y * invLen);
}

float Vec2f::cross(const Vec2f& right) const
{
    // Effectively the same as taking the determinant of a 2x2 matrix
    return x * right.y - y * right.x;
}

// ************************************************************

// ************************** Vec3f ***************************

Vec3f::Vec3f() : x(0), y(0), z(0) {}

Vec3f::Vec3f(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}

Vec3f::Vec3f(const Vec3f& other) : x(other.x), y(other.y), z(other.z) {}

Vec3f Vec3f::operator + (const Vec3f& right) const
{
    return Vec3f(x + right.x, y + right.y, z + right.z);
}

Vec3f Vec3f::operator - (const Vec3f& right) const
{
    return Vec3f(x - right.x, y - right.y, z - right.z);
}

Vec3f Vec3f::operator * (const float& scalar) const
{
    return Vec3f(x * scalar, y * scalar, z * scalar);
}

Vec3f Vec3f::operator / (const float& scalar) const
{
    return Vec3f(x / scalar, y / scalar, z / scalar);
}

float Vec3f::length() const
{
    return sqrtf(x * x + y * y + z * z);
}

float Vec3f::lengthSq() const
{
    return x * x + y * y + z * z;
}

float Vec3f::dot(const Vec3f& right) const
{
    return x * right.x + y * right.y + z * right.z;
}

Vec3f Vec3f::normalize() const
{
    float invLen = 1 / length();

    return Vec3f(x * invLen, y * invLen, z * invLen);
}

Vec3f Vec3f::cross(const Vec3f& right) const
{
    return Vec3f(
        y * right.z - z * right.y,
        z * right.x - x * right.z,
        x * right.y - y * right.x
    );
}

// ************************************************************

// ************************** Matrix44f ***************************

Matrix44f::Matrix44f() {}

Matrix44f::Matrix44f(std::initializer_list<float> values)
{
    if (values.size() != 16)
        throw std::invalid_argument("A 4 by 4 matrix must be initialized with 16 values.");

    const float* v = values.begin();
    for (int i = 0; i < 16; i++)
    {
        _data[i / 4][i % 4] = *(v++);
    }
}

Matrix44f Matrix44f::identity()
{
    return Matrix44f{1, 0, 0, 0, 
                    0, 1, 0, 0, 
                    0, 0, 1, 0, 
                    0, 0, 0, 1};
}

Matrix44f Matrix44f::encode(const float& radians, const Vec3f& axis, const Vec3f& translation)
{
    // TODO: implement this method
    return Matrix44f::identity();
}

const float* Matrix44f::operator [] (const int& idx) const
{
    return _data[idx];
}

float* Matrix44f::operator [] (const int& idx)
{
    return _data[idx];
}

Vec3f Matrix44f::operator * (const Vec3f& v) const
{
    float nx = v.x * _data[0][0] + v.y * _data[0][1] + v.z * _data[0][2] + _data[0][3];
    float ny = v.x * _data[1][0] + v.y * _data[1][1] + v.z * _data[1][2] + _data[1][3];
    float nz = v.x * _data[2][0] + v.y * _data[2][1] + v.z * _data[2][2] + _data[2][3];
    float nw = v.x * _data[3][0] + v.y * _data[3][1] + v.z * _data[3][2] + _data[3][3];

    if (nw != 1 && nw != 0)
    {
        nx /= nw;
        ny /= nw;
        nz /= nw;
    }

    return Vec3f(nx, ny, nz);
}

Matrix44f Matrix44f::operator * (const Matrix44f& right) const
{
    Matrix44f m;

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

Matrix44f Matrix44f::transpose() const
{
    Matrix44f transpose = *this;

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            transpose[j][i] = _data[i][j];

    return transpose;
}

Matrix44f Matrix44f::inverse() const
{
    // Uses Gauss-Jordan elimination

    Matrix44f inv = identity();
    Matrix44f m = *this;

    for (int c = 0; c < 4; c++)
    {
        // Step 1: Make sure pivot coefficient is not 0, if so swap rows
        if (m[c][c] == 0)
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
            float scalar = m[r][c] / m[c][c];

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
        float scalar = 1.0f / m[r][r];

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
            float scalar = m[r][c];

            for (int c2 = 0; c2 < 4; c2++)
            {
                m[r][c2] -= scalar * m[c][c2];
                inv[r][c2] -= scalar * inv[c][c2];
            }
        }
    }

    return inv;
}

// ************************************************************

// ************************** Triangle2D ***************************

Triangle2D::Triangle2D() {}

Triangle2D::Triangle2D(const Vec2f& v1, const Vec2f& v2, const Vec2f& v3) 
{
    v[0] = v1;
    v[1] = v2;
    v[2] = v3;
}

float Triangle2D::areaDoubled() const
{
    // NOTE: THE TRIANGLE'S POINTS MUST BE SPECIFIED IN COUNTER-CLOCKWISE ORDER!

    return (v[1] - v[0]).cross(v[2] - v[1]);
}

AABB2D Triangle2D::boundingBox() const
{
    return AABB2D(
        Vec2f(std::min(std::min(v[0].x, v[1].x), v[2].x), std::min(std::min(v[0].y, v[1].y), v[2].y)),
        Vec2f(std::max(std::max(v[0].x, v[1].x), v[2].x), std::max(std::max(v[0].y, v[1].y), v[2].y))
    );
}

bool Triangle2D::contains(const Vec2f& p) const
{
    // NOTE: THE TRIANGLE'S POINTS MUST BE SPECIFIED IN COUNTER-CLOCKWISE ORDER!
    // Uses the concept of determinant signs indicating rotational relationship (Pineda edge function)

    for (int i = 0; i < 3; i++)
    {
        Vec2f edge = v[(i + 1) % 3] - v[i];
        Vec2f diff = p - v[i];

        if (edge.cross(diff) < 0) return false;
    }

    return true;
}

// *****************************************************************

// ************************** AABB2D ***************************

AABB2D::AABB2D() : min(Vec2f()), max(Vec2f()) {}

AABB2D::AABB2D(const Vec2f& _min, const Vec2f& _max) : min(_min), max(_max) {}

bool AABB2D::contains(const Vec2f& p) const
{   
    return p.x >= min.x && p.x <= max.x && p.y >= min.y && p.y <= max.y;
}

// ************************************************************