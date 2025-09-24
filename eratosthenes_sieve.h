#ifndef HLIDAC_PES_ERASTHOTHENES_SIEVE_H
#define HLIDAC_PES_ERASTHOTHENES_SIEVE_H

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

        uint64_t find_next_set_bit(uint64_t start_bit_idx) const;
        void reset_bit(uint64_t bit_idx);

        static std::tuple<uint64_t, uint64_t> index_mask(uint64_t bit_idx);
        static std::tuple<uint64_t, uint64_t> index_reset_mask(uint64_t bit_idx);
        static uint64_t step(uint64_t prime_bit_idx);
        static uint64_t prime_idx_squared(uint64_t prime_bit_idx);
        static uint64_t prime(uint64_t prime_bit_idx);

        static uint64_t sieve_size_bits(uint64_t primes_to);
        static uint64_t sieve_size_bytes_alloc(uint64_t sieve_bits);
        static uint64_t init_pop_thread_count();

        struct prime_seed_metainfo {
            uint64_t sieve_step;
            uint64_t next_idx;
            std::vector<uint64_t> reset_masks;
        };

        void seed_bit_range(const std::vector<prime_seed_metainfo>& primes_meta,
            uint64_t num_of_patterns, uint64_t start_bit, uint64_t end_bit);
        static std::vector<uint64_t> seeding_pattern(uint64_t prime_bit_idx, uint64_t last_bit_idx, uint64_t sieve_step);

        // Helper constants for bit manipulation.
        static constexpr uint64_t BITS_PER_BYTE = 8;
        static constexpr uint64_t UINT64_BYTES = sizeof(uint64_t);
        static constexpr uint64_t UINT64_BITS = BITS_PER_BYTE * sizeof(uint64_t);
        static_assert(UINT64_BITS == 64, "Unexpected number of bits in uint64_t!");
        static constexpr uint64_t UINT64_ONE_MASK = 0xff'ff'ff'ff'ff'ff'ff'ff;
        static constexpr uint64_t UINT64_BIT_MASK = 0x00'00'00'00'00'00'00'3f;
        static constexpr uint64_t UINT64_IDX_SHIFT = 6ul;

        // Max number of uin64_t masks per small prime in the seeding vector.
        static constexpr uint64_t MAX_MASKS_PER_PRIME = 128;
        static_assert(MAX_MASKS_PER_PRIME >= 61, "At least 61 masks are required for the biggest seeding prime.");

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
