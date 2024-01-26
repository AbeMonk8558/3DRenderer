#include <initializer_list>
#include <stdexcept>
#include <immintrin.h>
#include "linearMath.hpp"

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

    void getPointer(float* ptr) const
    {
        _mm256_store_ps(ptr, data);
    }

    friend float_m256 max_m256(const float_m256& left, const float_m256& right);
    friend float_m256 min_m256(const float_m256& left, const float_m256& right);

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

Vec3f_m256 vec3fToVec3f_m256(const Vec3f* vecs)
{
    Vec3f_m256 res;

    return res;
}