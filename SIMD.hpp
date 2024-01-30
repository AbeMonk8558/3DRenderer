#pragma once

#include <initializer_list>
#include "linearMath.hpp"

namespace simd
{
    class float_m256;


    using Vec2f_m256 = Vec2<float_m256>;
    using Vec3f_m256 = Vec3<float_m256>;
    using Matrix44f_m256 = Matrix44<float_m256>;


    class float_m256
    {
    public:
        float_m256();

        float_m256(float value);

        float_m256(std::initializer_list<float> values);

        float_m256(const __m256& ps);

        float_m256 operator - () const;

        float_m256 operator + (const float_m256& right) const;

        float_m256 operator - (const float_m256& right) const;

        float_m256 operator * (const float_m256& right) const;

        float_m256 operator / (const float_m256& right) const;

        void operator -= (const float_m256& right);

        void operator += (const float_m256& right);

        void operator *= (const float_m256& right);

        void operator /= (const float_m256& right);

        __m256 operator < (const float_m256& right);

        __m256 operator <= (const float_m256& right);

        __m256 operator > (const float_m256& right);

        __m256 operator >= (const float_m256& right);

        float operator [] (int idx) const;

        float& operator [] (int idx);

        void storeData(float* ptr) const;
        
        friend float_m256 max_m256(const float_m256& left, const float_m256& right);
        friend float_m256 min_m256(const float_m256& left, const float_m256& right);
        friend float_m256 abs_m256(const float_m256& f);
        friend float_m256 sqrt_m256(const float_m256& f);

    private:
        __m256 data;
    };


    float_m256 max_m256(const float_m256& left, const float_m256& right);

    float_m256 min_m256(const float_m256& left, const float_m256& right);

    float_m256 abs_m256(const float_m256& f);

    float_m256 sqrt_m256(const float_m256& f);

    Vec3f_m256 vec3fToVec3f_m256(const Vec3f* vecs);

    Vec3f_m256 vec3fToVec3f_m256(const Vec3f* vecs, int num);
}

template <>
simd::float_m256 simd::Vec2f_m256::length() const;

template <>
simd::Vec3f_m256 simd::Matrix44f_m256::operator * (const simd::Vec3f_m256& v) const;

template <>
simd::Matrix44f_m256 simd::Matrix44f_m256::inverse() const;