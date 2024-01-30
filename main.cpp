#include "rasterizer.hpp"


/***************************** TODO LIST ******************************
1. Implement anti-aliasing (described in scratchapixel)
2. Implement SIMD optimization for rasterizer
3. Make vertex-line stage more accurate.
4. Encase the entire main file in a class so variables can be shared amongst functions.
5. Optimize Matrix44<float_m256> inverse method to be completely vectorized.
6. Implement edge case in Matrix44<float_m256> vector multiplication where w iz 0.
**********************************************************************/

int main(int argc, char** argv)
{
    Rasterizer::start();

    return 0;
}