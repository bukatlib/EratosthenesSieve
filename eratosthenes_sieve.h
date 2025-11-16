#ifndef HLIDAC_PES_ERASTHOTHENES_SIEVE_H
#define HLIDAC_PES_ERASTHOTHENES_SIEVE_H

#include <array>
#include <string>
#include <tuple>
#include <vector>
#include <stdint.h>
#include "cmake_vars.h"

class Eratosthenes {
    public:
        Eratosthenes(uint64_t primes_to);
        Eratosthenes(const std::string& filename);
		Eratosthenes(const Eratosthenes&) = delete;
		Eratosthenes& operator=(const Eratosthenes&) = delete;
		~Eratosthenes();

        uint64_t bit_size() const;
        uint64_t bytes_allocated() const;
        uint64_t to_sieved() const;

        void sieve_primes();

        uint64_t prime_count() const;
        bool is_prime(uint64_t prime_candidate) const;
        std::vector<uint64_t> primes() const;

        void write_primes_to_file(const std::string& filename) const;
        void write_image(const std::string& filename) const;

    private:
        // Max number of uin64_t masks per small prime in the seeding vector.
        static constexpr uint64_t MAX_MASKS_PER_PRIME = 128ul;
        static_assert(MAX_MASKS_PER_PRIME >= 61ul, "At least 61 masks are required for the biggest seeding prime.");

        // The maximal small prime number used for vectorized zeroing by using bit masking.
        // up to 41: 9488 bytes; up to 61: 13696 bytes (43, 47, 53, 59, 61 contain zero holes for 2x3x5)
        static constexpr uint64_t MAX_SMALL_PRIME = 59ul;

        // Constants for parallel processing.
        static constexpr uint64_t DEFAULT_SEGMENT_SIZE = 32ul * 1024ul;
        static constexpr uint64_t PARALLEL_SIEVE_MIN_BITS = 200e6;
        #if defined(__aarch64__)
        static constexpr uint64_t INIT_THREADS = 2ul;
        static constexpr uint64_t WORKCHUNK_SEGMENTS = 28ul;
        #else
        static constexpr uint64_t INIT_THREADS = 4ul;
        static constexpr uint64_t WORKCHUNK_SEGMENTS = 16ul;
        #endif

        // Core methods to get bit array idx and mask to extract or reset a bit.
        static std::tuple<uint64_t, uint64_t> index_mask(uint64_t bit_idx);
        static std::tuple<uint64_t, uint64_t> index_reset_mask(uint64_t bit_idx);
        void reset_bit(uint64_t bit_idx);

        // Helper bit methods to find the first set bit or check whether a number is in the wheel.
        uint64_t find_next_set_bit(uint64_t start_bit_idx) const;
        static bool is_in_image(uint64_t number, uint64_t* bit_idx = nullptr);
        static uint64_t sieve_bit_to_number(uint64_t bit_idx);

        // Calculation of the sieve size in bits and bytes.
        static uint64_t sieve_size_bits(uint64_t primes_to);
        static uint64_t sieve_size_elems(uint64_t sieve_bits);

        #ifdef WHEEL_2_3_5
        #pragma pack(push)
        #pragma pack(1)
        struct alignas(32) prime_wheel_steps {
            uint64_t step_idx : 3;  // Eight steps for 2x3x5 wheel.
            uint64_t bit_idx: 53;   // Will work up to 1000 TB of RAM
            uint32_t step0: 25;     // Will work up to 620 GB of RAM, or primes up to 20 * 10^12
            uint32_t step1: 25;
            uint32_t step2: 25;
            uint32_t step3: 25;
            uint32_t step4: 25;
            uint32_t step5: 25;
            uint32_t step6: 25;
            uint32_t step7: 25;
        };
        #pragma pack(pop)

        static constexpr std::array<uint64_t, WHEEL_STEPS> compute_wheel_2_3_5();
        static constexpr std::array<int8_t, WHEEL_CIRCUMFERENCE> compute_modulo_to_wheel_idx();
        #endif

        #ifdef WHEEL_2_3
        struct alignas(16) prime_wheel_steps {
            uint64_t step_idx: 1;  // Two steps for 2x3 wheel.
            uint64_t bit_idx: 63;
            uint32_t step0;
            uint32_t step1;
        };
        #endif

        #ifdef WHEEL_2
        struct prime_wheel_steps {
            uint64_t bit_idx;
            uint32_t step0;
        };
        #endif

        using prime_bit_masks = std::vector<std::vector<uint64_t>>;
        using prime_wheel_data = std::vector<prime_wheel_steps>;

        // List of numbers coprime to the wheel, and modulo to the wheel index (-1 if not in the wheel).
        alignas(64) static const std::array<uint64_t, WHEEL_STEPS> wheel;
        alignas(32) static const std::array<int8_t, WHEEL_CIRCUMFERENCE> modulo_to_idx;

        // Calculate the start bit index and steps for the given prime.
        static prime_wheel_steps wheel_steps(uint64_t prime);

        // Generate a reset masks for small primes (up to 61) used for vectorized sieving.
        static std::vector<uint64_t> seeding_pattern(uint64_t prime);

        // Sieving methods to reset bits for small, medium, and large steps.
        void sieve_bit_range(const prime_bit_masks& masks, const prime_wheel_data& wheel_data, uint64_t start_bit, uint64_t end_bit);
        void sieve_bit_range_small(const prime_bit_masks& masks, std::vector<uint64_t>& pattern_idxs, uint64_t start_bit, uint64_t end_bit);
        void sieve_bit_range_medium(prime_wheel_steps& prime_steps, uint64_t end_bit);
        void sieve_bit_range_large(prime_wheel_steps& prime_steps, uint64_t end_bit);

        // Exact number of bits required, number of allocated uint64_t elements, and the max prime to sieve.
        uint64_t sieve_bits;
        uint64_t sieve_size;
        uint64_t primes_to;

        // The sieve vector, odd numbers each consumes one bit.
		uint64_t *sieve;

        // Segment size, L1 data cache on ARM64, L2 data cache per hw thread on AMD64.
        uint64_t segment_size;
};

#endif
