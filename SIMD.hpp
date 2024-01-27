#pragma once

#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <immintrin.h>
#include "3DRenderer.hpp"

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

        bool operator == (const float_m256& right)
        {
            // For each element, sets every bit to 0 if not equal, otherwise to 1
            __m256 cmp = _mm256_cmp_ps(data, right.data, _CMP_EQ_OQ);

            // Movemask takes MSB for each float and stores it in a byte. If all are 1, the vectors are equal.
            return _mm256_movemask_ps(cmp) == 0xff;
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

    template <typename T>
    T absf(const T& val)
    {
        return fabsf(val);
    }

    template <>
    float_m256 absf<float_m256>(const float_m256& fv)
    {
        return abs_m256(fv);
    }

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
}