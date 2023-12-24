#pragma once

#include <initializer_list>

class Vec3f;
class Matrix44f;

class Vec3f
{
public:
    float x, y, z;

    Vec3f();
    Vec3f(float _x, float _y, float _z);
    Vec3f(const Vec3f& other);

    Vec3f operator + (const Vec3f& right) const;
    Vec3f operator - (const Vec3f& right) const;
    Vec3f operator * (const float& scalar) const;
    Vec3f operator / (const float& scalar) const;

    float length() const;
    float lengthSq() const;
    float dot(const Vec3f& right) const;
    Vec3f normalize() const;
    Vec3f cross(const Vec3f& right) const;
};

class Matrix44f
{
public:
    Matrix44f();
    Matrix44f(std::initializer_list<float> values);
    static Matrix44f identity();
    static Matrix44f encode(const float& radians, const Vec3f& axis, const Vec3f& translation);

    const float* operator [] (const int& idx) const;
    float* operator [] (const int& idx);
    Vec3f operator * (const Vec3f& v) const;
    Matrix44f operator * (const Matrix44f& right) const;

    Matrix44f transpose() const;
    Matrix44f inverse() const;

private:
    float _data[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
};