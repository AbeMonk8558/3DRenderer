#pragma once

#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <immintrin.h>
#include "3DRenderer.hpp"
#include "linearMath.hpp"

namespace simd
{
    class float_m256
    {
    public:
        float_m256() : data(_mm256_set1_ps(0.0f)) {}

        float_m256(float value) : data(_mm256_set1_ps(value)) {}

        float_m256(std::initializer_list<float> values) 
        {
            if (values.size() != 8)
                throw std::invalid_argument("A 256-bit SIMD vector must be initialized with 8x 4-byte values.");

            data = _mm256_load_ps(values.begin());
        }

        float_m256(const __m256& ps) : data(ps) {}

        float_m256 operator - () const
        {
            // Represents a float with every bit 0 except the MSB, which is 1
            int i = 0x80000000;
            float* fMaskVal = reinterpret_cast<float*>(&i);

            __m256 xorMask = _mm256_set1_ps(*fMaskVal);
            return _mm256_xor_ps(data, xorMask);
        }

        float_m256 operator + (const float_m256& right) const
        {
            return float_m256(_mm256_add_ps(data, right.data));
        }

        float_m256 operator - (const float_m256& right) const
        {
            return float_m256(_mm256_sub_ps(data, right.data));
        }

        float_m256 operator * (const float_m256& right) const
        {
            return float_m256(_mm256_mul_ps(data, right.data));
        }

        float_m256 operator / (const float_m256& right) const
        {
            return float_m256(_mm256_div_ps(data, right.data));
        }

        void operator -= (const float_m256& right)
        {
            data = _mm256_sub_ps(data, right.data);
        }

        void operator += (const float_m256& right)
        {
            data = _mm256_add_ps(data, right.data);
        }

        void operator *= (const float_m256& right)
        {
            data = _mm256_mul_ps(data, right.data);
        }

        void operator /= (const float_m256& right)
        {
            data = _mm256_div_ps(data, right.data);
        }

        __m256 operator < (const float_m256& right)
        {
            return _mm256_cmp_ps(data, right.data, _CMP_LT_OQ);
        }

        __m256 operator <= (const float_m256& right)
        {
            return _mm256_cmp_ps(data, right.data, _CMP_LE_OQ);
        }

        __m256 operator > (const float_m256& right)
        {
            return _mm256_cmp_ps(data, right.data, _CMP_GT_OQ);
        }

        __m256 operator >= (const float_m256& right)
        {
            return _mm256_cmp_ps(data, right.data, _CMP_GE_OQ);
        }

        float operator [] (int idx) const
        {
            return data[idx];
        }

        float& operator [] (int idx)
        {
            return data[idx];
        }

        void storeData(float* ptr) const
        {
            _mm256_store_ps(ptr, data);
        }
        
        friend float_m256 max_m256(const float_m256& left, const float_m256& right);
        friend float_m256 min_m256(const float_m256& left, const float_m256& right);
        friend float_m256 abs_m256(const float_m256& f);
        friend float_m256 sqrt_m256(const float_m256& f);

    private:
        __m256 data;
    };

    float_m256 max_m256(const float_m256& left, const float_m256& right)
    {
        return float_m256(_mm256_max_ps(left.data, right.data));
    }

    float_m256 min_m256(const float_m256& left, const float_m256& right)
    {
        return float_m256(_mm256_min_ps(left.data, right.data));
    }

    float_m256 abs_m256(const float_m256& f)
    {
        // Represents a float with every bit 1 except the MSB, which is 0
        int i = 0x7FFFFFFF;
        float* fMaskVal = reinterpret_cast<float*>(&i);

        __m256 andMask = _mm256_set1_ps(*fMaskVal);
        return _mm256_and_ps(f.data, andMask);
    }

    float_m256 sqrt_m256(const float_m256& f)
    {
        return float_m256(_mm256_sqrt_ps(f.data));
    }

    // template <typename T>
    // T absf(const T& val)
    // {
    //     return fabsf(val);
    // }

    // template <>
    // float_m256 absf<float_m256>(const float_m256& fv)
    // {
    //     return abs_m256(fv);
    // }

    Vec3f_m256 vec3fToVec3f_m256(const Vec3f* vecs)
    {
        Vec3f_m256 res;
        res.x = {vecs[0].x, vecs[1].x, vecs[2].x, vecs[3].x, vecs[4].x, vecs[5].x, vecs[6].x, vecs[7].x};
        res.y = {vecs[0].y, vecs[1].y, vecs[2].y, vecs[3].y, vecs[4].y, vecs[5].y, vecs[6].y, vecs[7].y};
        res.z = {vecs[0].z, vecs[1].z, vecs[2].z, vecs[3].z, vecs[4].z, vecs[5].z, vecs[6].z, vecs[7].z};

        return res;
    }

    Vec3f_m256 vec3fToVec3f_m256(const Vec3f* vecs, int num)
    {
        Vec3f_m256 res;

        // A value of 0 would lead to zero-division errors
        res.x = _mm256_set1_ps(1.0f);
        res.y = _mm256_set1_ps(1.0f);
        res.z = _mm256_set1_ps(1.0f);

        for (int i = 0; i < num; i++)
        {
            res.x[i] = vecs[i].x;
            res.y[i] = vecs[i].y;
            res.z[i] = vecs[i].z;
        }

        return res;
    }

    Vec2f_m256 vec3fToUniformVec3f_m256(const Vec3f& vec)
    {
        
    }
}

template <>
simd::float_m256 simd::Vec2f_m256::length() const
{
    return sqrt_m256(x * x + y * y);
}

template <>
simd::Vec3f_m256 simd::Matrix44f_m256::operator * (const simd::Vec3f_m256& v) const
{
    simd::float_m256 nx = v.x * _data[0][0] + v.y * _data[0][1] + v.z * _data[0][2] + _data[0][3];
    simd::float_m256 ny = v.x * _data[1][0] + v.y * _data[1][1] + v.z * _data[1][2] + _data[1][3];
    simd::float_m256 nz = v.x * _data[2][0] + v.y * _data[2][1] + v.z * _data[2][2] + _data[2][3];
    simd::float_m256 nw = v.x * _data[3][0] + v.y * _data[3][1] + v.z * _data[3][2] + _data[3][3];

    return simd::Vec3f_m256(nx, ny, nz);
}

template <>
simd::Matrix44f_m256 simd::Matrix44f_m256::inverse() const
{
    simd::Matrix44f_m256 inv = identity();
    simd::Matrix44f_m256 m = *this;

    for (int c = 0; c < 4; c++)
    {
        for (int i = 0; i < 8; i++)
        {
            if (m[c][c][i] == 0)
            {
                int max = c; // Pivot coefficient has row index equal to column index

                // Find row with maximum absolute value coefficient in same column
                for (int r = 0; r < 4; r++)
                    if (fabsf(m[r][c][i]) > fabsf(m[max][c][i])) max = r;

                if (max == c) return identity(); // TODO: Should probably throw exception or something

                for (int c2 = 0; c2 < 4; c2++)
                {
                    std::swap(m[max][c2][i], m[c][c2][i]);
                    std::swap(inv[max][c2][i], inv[c][c2][i]);
                }
            }
        }

        // Step 2: Perform forward substitution on the column, setting all coefficients below
        //         the pivot within the column to 0 through row addition and scaling
        for (int r = c + 1; r < 4; r++) // Start at the row below the pivot in the column
        {
            simd::float_m256 scalar = m[r][c] / m[c][c];

            for (int c2 = 0; c2 < 4; c2++)
            {
                m[r][c2] -= scalar * m[c][c2];
                inv[r][c2] -= scalar * inv[c][c2];
            }

            m[r][c] = simd::float_m256(0); // Just to be safe
        }
    }

    // Step 3: Set all the pivot coefficients to 1 through row scaling
    for (int r = 0; r < 4; r++)
    {
        simd::float_m256 scalar = simd::float_m256(1) / m[r][r];

        for (int c = 0; c < 4; c++)
        {
            m[r][c] *= scalar;
            inv[r][c] *= scalar;
        }

        m[r][r] = simd::float_m256(1); // Just to be safe
    }

    // Step 4: Perform backward substitution on all columns above pivot
    for (int c = 0; c < 4; c++)
    {
        for (int r = 0; r < c; r++)
        {
            simd::float_m256 scalar = m[r][c];

            for (int c2 = 0; c2 < 4; c2++)
            {
                m[r][c2] -= scalar * m[c][c2];
                inv[r][c2] -= scalar * inv[c][c2];
            }
        }
    }

    return inv;
}