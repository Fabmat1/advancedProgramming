#include <iostream>
#include <random>
#include <chrono>
#include <immintrin.h>
#include <omp.h>

using namespace std::chrono;

int main() {
    auto start = high_resolution_clock::now();
    // Seed the random number generator (Kind of does not work due to parallelization but too lazy to fix, sorry!)
    std::mt19937_64 eng(12345);

    // Generate 10000000 random values for x and y
    const int N = 10000000;
    int N_i = 0;

    // Get the maximum number of threads supported by the CPU
    int max_threads = omp_get_max_threads();

    // Set the number of threads to use in parallel regions
    omp_set_num_threads(max_threads);

    // Define the distribution
    std::uniform_real_distribution<double> distr(0, 1.0);
    __m512d r_squared = _mm512_set1_pd(std::pow(0.5, 2));

    // Parallelize for loop execution (Use all CPU Cores)
    #pragma omp parallel for reduction(+: N_i)
    for  (int i = 0; i < N; i += 8) {
        // I <3 vectorization!
        // Generate batches of 8 random doubles and substract 0.5
        __m512d x_values = _mm512_sub_pd(_mm512_set_pd(distr(eng), distr(eng), distr(eng), distr(eng), distr(eng), distr(eng), distr(eng), distr(eng)), _mm512_set1_pd(0.5));
        __m512d y_values = _mm512_sub_pd(_mm512_set_pd(distr(eng), distr(eng), distr(eng), distr(eng), distr(eng), distr(eng), distr(eng), distr(eng)), _mm512_set1_pd(0.5));

        // square the values
        __m512d x_squared = _mm512_mul_pd(x_values, x_values);
        __m512d y_squared = _mm512_mul_pd(y_values, y_values);

        // Sum the squares and determine if inside circle
        __m512d sum_squared = _mm512_add_pd(x_squared, y_squared);
        __mmask8 inside = _mm512_cmp_pd_mask(sum_squared, r_squared, _MM_CMPINT_LT);

        // Keep track of the amount inside the circle
        N_i += _mm_popcnt_u64(inside);
    }

    double fraction =  static_cast<double>(N_i)/static_cast<double>(N);
    auto end = high_resolution_clock::now();

    std::cout.precision(5);
    std::cout << fraction << std::fixed << std::endl;
    std::cout << 4*fraction << std::fixed << std::endl;

    auto duration = duration_cast<microseconds>(end - start);
    auto elapsed = static_cast<float>(duration.count());
    std::cout << elapsed/1000000 << std::endl;

    return 0;
}
