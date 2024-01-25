#pragma once

#include <initializer_list>
#include <immintrin.h>
#include <raylib.h>

class Vec2f;
class Vec3f;
class Matrix44f;
// class Triangle2D;
// class AABB2D;

class Vec2f
{
public:
    float x, y;

    Vec2f();
    Vec2f(float _x, float _y);
    Vec2f(const Vec2f& other);
    Vec2f(const Vector2& raylibVec);

    Vec2f operator + (const Vec2f& right) const;
    Vec2f operator - (const Vec2f& right) const;
    Vec2f operator * (const float& scalar) const;
    Vec2f operator / (const float& scalar) const;

    float length() const;
    float lengthSq() const;
    float dot(const Vec2f& right) const;
    Vec2f normalize() const;
    float cross(const Vec2f& right) const;
};

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

// class Triangle2D
// {
// public:
//     Vec2f v[3] = {Vec2f(), Vec2f(), Vec2f()}; // Must be in counter-clockwise order

//     Triangle2D();
//     Triangle2D(const Vec2f& v1, const Vec2f& v2, const Vec2f& v3);

//     float areaDoubled() const;
//     AABB2D boundingBox() const;
//     bool contains(const Vec2f& p) const;
// };

// class AABB2D // Axis-aligned bounding box
// {
// public:
//     Vec2f min, max;

//     AABB2D();
//     AABB2D(const Vec2f& _min, const Vec2f& _max);

//     bool contains(const Vec2f& p) const;
// };

