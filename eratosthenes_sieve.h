#ifndef HLIDAC_PES_ERASTHOTHENES_SIEVE_H
#define HLIDAC_PES_ERASTHOTHENES_SIEVE_H

#include <array>
#include <string>
#include <tuple>
#include <vector>
#include <stdint.h>

class Eratosthenes {
    public:
        Eratosthenes(uint64_t primes_to);
		Eratosthenes(const Eratosthenes&) = delete;
		Eratosthenes& operator=(const Eratosthenes&) = delete;
		~Eratosthenes();

        uint64_t bit_size() const;
        uint64_t bytes_allocated() const;

        void sieve_primes();
        uint64_t prime_count();
        bool is_prime(uint64_t prime_candidate) const;
        std::vector<uint64_t> primes() const;
        void write_primes_to_file(const std::string& filename) const;
        void write_image(const std::string& filename) const;

    private:
        // Helper constants for bit manipulation.
        static constexpr uint64_t BITS_PER_BYTE = 8ul;
        static constexpr uint64_t UINT64_BYTES = sizeof(uint64_t);
        static constexpr uint64_t UINT64_BITS = BITS_PER_BYTE * sizeof(uint64_t);
        static_assert(UINT64_BITS == 64ul, "Unexpected number of bits in uint64_t!");
        static constexpr uint64_t UINT64_ONE_MASK = 0xff'ff'ff'ff'ff'ff'ff'ff;
        static constexpr uint64_t UINT64_BIT_MASK = 0x00'00'00'00'00'00'00'3f;
        static constexpr uint64_t UINT64_IDX_SHIFT = 6ul;

        // Max number of uin64_t masks per small prime in the seeding vector.
        static constexpr uint64_t MAX_MASKS_PER_PRIME = 128ul;
        static_assert(MAX_MASKS_PER_PRIME >= 61ul, "At least 61 masks are required for the biggest seeding prime.");

        // The maximal small prime number for vectorized masking (43, 47, 53, 59, 61 contain zero holes in the pattern).
        // up to 41: 9488 bytes; up to 61: 13696 bytes
        static constexpr uint64_t MAX_SMALL_PRIME = 41ul;

        // Wheel 2x3x5 (basis) helper constants, computation functions, read-only wheel data.
        static constexpr uint64_t WHEEL_NUMBER_OF_COPRIMES = 8ul;
        static constexpr uint64_t WHEEL_BITS = 8ul;
        static constexpr uint64_t WHEEL_CIRCUMFERENCE = 30ul;

        static constexpr std::array<uint64_t, WHEEL_NUMBER_OF_COPRIMES> compute_wheel_2_3_5();
        static constexpr std::array<int8_t, WHEEL_CIRCUMFERENCE> compute_modulo_to_wheel_idx();

        // List of numbers coprime to wheel 2x3x5 (basis), align to a typical cache line.
        alignas(64) static const std::array<uint64_t, WHEEL_NUMBER_OF_COPRIMES> wheel;

        // Mapping from a modulo 2x3x5 to the wheel index (-1 if not in the wheel).
        alignas(32) static const std::array<int8_t, WHEEL_CIRCUMFERENCE> modulo_to_idx;


        uint64_t find_next_set_bit(uint64_t start_bit_idx) const;
        void reset_bit(uint64_t bit_idx);

        static std::tuple<uint64_t, uint64_t> index_mask(uint64_t bit_idx);
        static std::tuple<uint64_t, uint64_t> index_reset_mask(uint64_t bit_idx);
        static uint64_t sieve_bit_to_number(uint64_t bit_idx);

        static bool is_in_image(uint64_t number, uint64_t* bit_idx = nullptr);
        static std::array<uint64_t, 8> wheel_modulo_init_idxs(uint64_t prime);

        static uint64_t sieve_size_bits(uint64_t primes_to);
        static uint64_t sieve_size_bytes_alloc(uint64_t sieve_bits);
        static uint64_t init_pop_thread_count();

        struct prime_seed_metainfo {
            uint64_t prime;
            std::array<uint64_t, 8> modulo_init_idxs;
            std::vector<uint64_t> reset_masks;
        };

        void seed_bit_range(const std::vector<prime_seed_metainfo>& primes_meta,
            uint64_t num_of_patterns, uint64_t start_bit, uint64_t end_bit);
        static std::vector<uint64_t> seeding_pattern(uint64_t prime);

        // Exact number of bits required, allocated memory in bytes (rounded to 8-byte multiple). and the sieve vector.
        uint64_t sieve_bits;
        uint64_t sieve_bytes_alloc;

        // Number of found prime numbers and the value to which primes were found.
        uint64_t found_primes;
        uint64_t primes_to;

        // The sieve vector, odd numbers each consumes one bit.
		uint64_t *sieve;
		uint64_t sieve_size;

        // Segment size, L1 data cache on ARM64, L2 data cache per hw thread on AMD64.
        uint64_t segment_size;
        static constexpr uint64_t WORKCHUNK_SEGMENTS = 16ul;
};

#endif
