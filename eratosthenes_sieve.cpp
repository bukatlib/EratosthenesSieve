#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "cpu_info.h"
#include "eratosthenes_sieve.h"

using namespace std;
using namespace std::chrono;

Eratosthenes::Eratosthenes(uint64_t primes_to) : 
    sieve_bits(sieve_size_bits(primes_to)), sieve_bytes_alloc(sieve_size_bytes_alloc(sieve_bits)),
    found_primes(0ul), primes_to(primes_to), sieve_size(sieve_bytes_alloc >> 3)
{
	sieve = new uint64_t[sieve_size];

    uint64_t num_of_threads = Eratosthenes::init_pop_thread_count();
    #pragma omp parallel for if(sieve_size > 10000000) num_threads(num_of_threads)
	for (uint64_t i = 0; i < sieve_size; ++i)
		sieve[i] = UINT64_ONE_MASK;

    // Set the unused tail bits to false to avoid primes behind the requested range.
    uint64_t last_uint64_valid_bits = sieve_bits % UINT64_BITS;
    if (last_uint64_valid_bits > 0) {
        uint64_t reset_mask = (1ul << last_uint64_valid_bits) - 1ul;
        sieve[sieve_bits>>UINT64_IDX_SHIFT] &= reset_mask;
    }

    CpuInfo cpuinfo;
    segment_size = 32ul * 1024ul;
    prop_ca architecture = cpuinfo.arch();
    if (architecture.has_value())   {
        prop_ui cache_size;
        switch (architecture.value())   {
            case CpuArch::AArch64:
                // Use the exclusive L1 data cache size (best on RPI5).
                cache_size = cpuinfo.cache_size_exclusive(1, CacheType::Data);
                if (cache_size)
                    segment_size = cache_size.value();
                break;
            case CpuArch::AMD64:
                // Use 75% of the exclusive L2 data cache size (best on Ryzen 7950x).
                cache_size = cpuinfo.cache_size_exclusive(2, CacheType::Data);
                if (cache_size)	{
					segment_size = cache_size.value();
                    segment_size -= cache_size.value() >> 3ul;
				}
                break;
        }
    }

    // Segment size is in bits, adjust the value accordingly.
    segment_size *= BITS_PER_BYTE;

    // Prime squared index is up to 1859 for masks, avoid reseting the previous bits.
    assert(segment_size >= 2048ul);
    // Segment size should be a multiple of uint64_t size!
    assert(segment_size % UINT64_BITS == 0);
}

Eratosthenes::~Eratosthenes()	{
	delete[] sieve;
}

uint64_t Eratosthenes::bit_size() const {
    return sieve_bits;
}

uint64_t Eratosthenes::bytes_allocated() const {
    return sieve_bytes_alloc;
}

uint64_t Eratosthenes::init_pop_thread_count()   {
    #if defined(__aarch64__)
    return 2ul;
    #else
    return 4ul;
    #endif
}

void Eratosthenes::sieve_primes()    {
    // Handle the case when sieve is empty.
    if (sieve_bits == 0ul)
        return;

    uint64_t prime_bit_idx = 0;
    vector<prime_seed_metainfo> primes_meta;
    uint64_t segment_start = 0ul, segment_end = min(segment_size, sieve_bits);

    uint64_t patterns = 0ul;
    bool parallel_sieve = sieve_bits >= 200e6;

    while ((prime_bit_idx = find_next_set_bit(prime_bit_idx)) < sieve_bits)  {
        if (prime_bit_idx >= segment_end)   {
            segment_start = segment_end;
            segment_end = min(segment_start + segment_size, sieve_bits);
            seed_bit_range(primes_meta, patterns, segment_start, segment_end);
        }

        assert(prime_bit_idx < segment_end);

        uint64_t prime_squared_idx = prime_idx_squared(prime_bit_idx);
        uint64_t bit_idx = prime_squared_idx;
        uint64_t sieve_step = step(prime_bit_idx);

        if (prime_squared_idx < sieve_bits) {
            // Start from the square of the prime (lower had to be already seeded).
            while (bit_idx < segment_end)   {
                reset_bit(bit_idx);
                bit_idx += sieve_step;
            }
        } else {
            // Seed only with the primes up to the root square of n.
            break;
        }

        if (sieve_step < UINT64_BITS && bit_idx < sieve_bits)   {
            assert(prime_squared_idx < segment_size);
            primes_meta.push_back({sieve_step, 0ul, seeding_pattern(prime_bit_idx, bit_idx, sieve_step)});
            ++patterns;
        } else {
            primes_meta.push_back({sieve_step, bit_idx, {}});
        }
        
        prime_bit_idx++;
    }

    uint64_t chunk_size = WORKCHUNK_SEGMENTS * segment_size;
    uint64_t chunk_last_idx = segment_end + chunk_size * ((sieve_bits - segment_end) / chunk_size);

    #pragma omp parallel if(parallel_sieve)
    {
        #pragma omp for schedule(dynamic)
        for (uint64_t b = segment_end; b < chunk_last_idx; b += chunk_size)   {
            uint64_t start_bit = b, end_bit = b + chunk_size;
            seed_bit_range(primes_meta, patterns, start_bit, end_bit);
        }

        #pragma omp for schedule(static)
        for (uint64_t b = chunk_last_idx; b < sieve_bits; b += segment_size)    {
            uint64_t segment_start = b, segment_end = min(sieve_bits, b + segment_size);
            seed_bit_range(primes_meta, patterns, segment_start, segment_end);
        }
    }
}

uint64_t Eratosthenes::prime_count()  {
    if (sieve_size > 0ul)  {
        uint64_t size_aligned = sieve_size & 0xffffffff'fffffffc;
        uint64_t thread_count = Eratosthenes::init_pop_thread_count();
        uint64_t count0 = 0ul, count1 = 0ul, count2 = 0ul, count3 = 0ul, count4 = 0ul;

        #pragma omp parallel for reduction(+: count0,count1,count2,count3) if(sieve_size > 10000000) num_threads(thread_count)
        for (uint64_t uint64_idx = 0; uint64_idx < size_aligned; uint64_idx += 4)  {
            count0 += popcount(sieve[uint64_idx + 0]);
            count1 += popcount(sieve[uint64_idx + 1]);
            count2 += popcount(sieve[uint64_idx + 2]);
            count3 += popcount(sieve[uint64_idx + 3]);
        }

        for (uint64_t uint64_idx = size_aligned; uint64_idx < sieve_size; ++uint64_idx)
            count4 += popcount(sieve[uint64_idx]);

        // Add one for prime 2.
        found_primes = 1ul + count0 + count1 + count2 + count3 + count4;
    } else if (primes_to == 2)  {
        found_primes = 1ul;
    } else {
        found_primes = 0ul;
    }

    return found_primes;
}

vector<uint64_t> Eratosthenes::primes() const {
    vector<uint64_t> primes;
    primes.reserve(found_primes);
    if (primes_to >= 2ul)
        primes.push_back(2ul);
    
    if (sieve_size > 0ul)  {
        uint64_t prime_bit_idx = 0;
        while ((prime_bit_idx = find_next_set_bit(prime_bit_idx)) < sieve_bits) {
            primes.push_back(prime(prime_bit_idx));
            ++prime_bit_idx;
        }
    }

    return primes;
}

void Eratosthenes::write_primes_to_file(const string& filename) const {
    // Bigger C++ stream buffer to improve write performance.
    constexpr uint64_t buffer_size = 1ul<<20ul;
    char iobuffer[buffer_size];

    ofstream OUT;
    OUT.rdbuf()->pubsetbuf(iobuffer, buffer_size);
    OUT.open(filename.c_str(), ios::out | ios::trunc);
    if (!OUT.good()) {
    	cerr<<"Failed to create the file to write primes!"<<endl;
        return;
    }

    if (primes_to >= 2ul)
        OUT<<2<<"\n";

    if (sieve_size > 0ul)  {
        array<uint64_t, 4> prime_buffer;
        uint64_t buffer_idx = 0ul, prime_bit_idx = 0ul;
        while ((prime_bit_idx = find_next_set_bit(prime_bit_idx)) < sieve_bits) {
            uint64_t prime_val = prime(prime_bit_idx);

            if (buffer_idx == prime_buffer.size())   {
                buffer_idx = 0ul;
                string n0 = to_string(prime_buffer[0]), n1 = to_string(prime_buffer[1]);
                string n2 = to_string(prime_buffer[2]), n3 = to_string(prime_buffer[3]);
                OUT<<n0<<"\n"<<n1<<"\n"<<n2<<"\n"<<n3<<"\n";
            }

            prime_buffer[buffer_idx++] = prime_val;
            
            ++prime_bit_idx;
        }

        for (uint64_t i = 0ul; i < buffer_idx; ++i)
            OUT<<prime_buffer[i]<<"\n";
    }

    OUT.flush();
    OUT.close();
}

void Eratosthenes::write_image(const string& filename) const {
    ofstream OUT(filename.c_str(), ios::out | ios::trunc | ios::binary);

    for (uint64_t i = 0ul; i < sieve_size; ++i)
        OUT.write(reinterpret_cast<char*>(&sieve[i]), UINT64_BYTES);

    OUT.flush();
    OUT.close();
}

bool Eratosthenes::is_prime(uint64_t prime_candidate) const {
    if (prime_candidate <= 1ul)  {
        return false;
    } else if (prime_candidate == 2ul)  {
        return true;
    } else if ((prime_candidate & 0x01) == 0ul) {
        return false;
    } else if (prime_candidate <= primes_to)    {
        uint64_t bit_idx = (prime_candidate - 3ul) >> 1ul;
        auto [idx, mask] = index_mask(bit_idx);
        return (sieve[idx] & mask) > 0ul;
    } else {
        throw out_of_range("Increase the size of the sieve, too large prime to check!");
    }
}

uint64_t Eratosthenes::find_next_set_bit(uint64_t start_bit_idx) const {
    uint64_t uint64_idx = start_bit_idx >> UINT64_IDX_SHIFT;
    uint64_t first_bit_idx = start_bit_idx - UINT64_BITS * uint64_idx;
    assert(first_bit_idx < UINT64_BITS);
    uint64_t mask_first = UINT64_ONE_MASK<<first_bit_idx;

    assert(start_bit_idx <= sieve_bits);

    uint64_t right_zero_bits = 0ul;
    uint64_t uint64_val = sieve[uint64_idx] & mask_first;
    
    while ((right_zero_bits = countr_zero(uint64_val)) == UINT64_BITS && ++uint64_idx < sieve_size)  {
        uint64_val = sieve[uint64_idx];
    }

    if (uint64_idx < sieve_size)
        // Number of right bit zeros is lower than 64, we found a set bit in the current uint64_t.
        return UINT64_BITS * uint64_idx + right_zero_bits;
    else
        // Return the size to indicate the end of sieve (not found indicator).
        return sieve_bits;
}

void Eratosthenes::reset_bit(uint64_t bit_idx) {
    // Reset the bit in the sieve.
    assert(bit_idx < sieve_bits);
    #if defined(__aarch64__)
    auto [uint64_idx, bit_mask] = index_mask(bit_idx);
    sieve[uint64_idx] &= ~bit_mask;
    #else
    auto [uint64_idx, reset_bit_mask] = index_reset_mask(bit_idx);
    sieve[uint64_idx] &= reset_bit_mask;
    #endif
}

tuple<uint64_t, uint64_t> Eratosthenes::index_mask(uint64_t bit_idx)   {
    // Index to sieve vector and set/get bitmask.
    return {bit_idx >> UINT64_IDX_SHIFT, 1ul<<(UINT64_BIT_MASK & bit_idx)};
}

tuple<uint64_t, uint64_t> Eratosthenes::index_reset_mask(uint64_t bit_idx)   {
    // Index to sieve vector and set/get bitmask.
    return {bit_idx >> UINT64_IDX_SHIFT, rotl(0xffff'ffff'ffff'fffe, UINT64_BIT_MASK & bit_idx)};
}

uint64_t Eratosthenes::step(uint64_t prime_bit_idx)  {
    // Calculate the seeding step for the prime.
    return 2ul * prime_bit_idx + 3ul;
}

uint64_t Eratosthenes::prime_idx_squared(uint64_t prime_bit_idx)   {
    // If prime_bit_idx points to a prime, then the return idx points to the power of the prime.
    return 2 * prime_bit_idx * prime_bit_idx + 6 * prime_bit_idx + 3;
}

uint64_t Eratosthenes::prime(uint64_t prime_bit_idx)   {
    // Calculate the prime that corresponds to the bit index, the same as step value.
    return 2ul * prime_bit_idx + 3ul;
}

uint64_t Eratosthenes::sieve_size_bits(uint64_t primes_to)  {
    // Only odd numbers are in the sieve from 3 upwards.
    if (primes_to >= 3)
        return ((primes_to - 3) >> 1) + 1;
    else
        return 0;
}

uint64_t Eratosthenes::sieve_size_bytes_alloc(uint64_t sieve_bits)  {
    // Round up the size to 64 bit multiple.
    return UINT64_BYTES * ((sieve_bits + UINT64_BITS - 1) / UINT64_BITS);
}

void Eratosthenes::seed_bit_range(const vector<prime_seed_metainfo>& primes_meta,
                    uint64_t num_of_patterns, uint64_t start_bit, uint64_t end_bit)    {

    vector<uint64_t> next_idx(primes_meta.size(), 0ul);

    for (uint64_t i = 0ul; i < num_of_patterns; ++i)  {
        assert(start_bit >= segment_size);
        uint64_t mask_count = primes_meta[i].reset_masks.size();
        assert(mask_count > 0ul);
        uint64_t pattern_idx = ((start_bit - segment_size) / UINT64_BITS) % mask_count;
        next_idx[i] = pattern_idx;
    }

    bool terminate = false;
    uint64_t prime_idx = num_of_patterns;

    while (prime_idx < primes_meta.size() && !terminate)    {
        uint64_t init_bit_idx = primes_meta[prime_idx].next_idx;
        uint64_t sieve_step = primes_meta[prime_idx].sieve_step;
        uint64_t min_prime_bit = max(start_bit, init_bit_idx);
        uint64_t cur_bit_idx = init_bit_idx + ((min_prime_bit - init_bit_idx + sieve_step - 1ul) / sieve_step) * sieve_step;
        next_idx[prime_idx] = cur_bit_idx;
        terminate = init_bit_idx >= end_bit;
        ++prime_idx;
    }

    for (uint64_t segment_start = start_bit; segment_start < end_bit; segment_start += segment_size)    {
        uint64_t segment_end = min(end_bit, segment_start + segment_size); 

        // Use pattern to seed small primes.
        for (uint64_t i = 0ul; i < num_of_patterns; ++i) {
            const vector<uint64_t>& pattern = primes_meta[i].reset_masks;
            uint64_t pattern_chunks = pattern.size();
            uint64_t segment_bits = segment_end - segment_start;
            uint64_t partial_bits = segment_bits & UINT64_BIT_MASK;
            uint64_t sieve_idx = (segment_start >> UINT64_IDX_SHIFT), pattern_idx = next_idx[i];
            uint64_t pattern_writes = (segment_bits >> UINT64_IDX_SHIFT) + (partial_bits == 0ul ? 0ul : 1ul) ;

            uint64_t partial_pattern_count = min(pattern_writes, pattern_chunks - pattern_idx);
            for (uint64_t i = 0; i < partial_pattern_count; ++i)
                sieve[sieve_idx++] &= pattern[pattern_idx + i];
            pattern_writes -= partial_pattern_count;

            while (pattern_writes >= pattern_chunks) {
                for (uint64_t i = 0; i < pattern_chunks; ++i)
                    sieve[sieve_idx++] &= pattern[i];
                pattern_writes -= pattern_chunks;
            }

            for (uint64_t i = 0; i < pattern_writes; ++i)
                sieve[sieve_idx++] &= pattern[i];

            next_idx[i] = pattern_writes;
        }

        // Reset idividual bits for large primes.
        for (uint64_t i = num_of_patterns; i < prime_idx; ++i) {
            uint64_t one_step = primes_meta[i].sieve_step;
            uint64_t two_steps = one_step << 1ul, four_steps = one_step << 2ul;
            uint64_t three_steps = two_steps + one_step;
            assert(segment_end >= four_steps);  // only 3^2 is less than 4*3 should not happen

            uint64_t bit_idx = next_idx[i], unroll_max_idx = segment_end - four_steps;

            while (bit_idx < unroll_max_idx)    {
                reset_bit(bit_idx);
                reset_bit(bit_idx + one_step);
                reset_bit(bit_idx + two_steps);
                reset_bit(bit_idx + three_steps);
                bit_idx += four_steps;
            }

            while (bit_idx < segment_end)  {
                reset_bit(bit_idx);
                bit_idx += one_step;
            }

            next_idx[i] = bit_idx;
        }
    }
}

vector<uint64_t> Eratosthenes::seeding_pattern(uint64_t prime_bit_idx, uint64_t last_bit_idx, uint64_t sieve_step) {
    // Calculate a seeding pattern for a small prime number.
    vector<uint64_t> masks;
    uint64_t prime_number = prime(prime_bit_idx);
    uint64_t repeats = MAX_MASKS_PER_PRIME / prime_number;
    assert(repeats >= 1ul);
    masks.reserve(prime_number*repeats);

    uint64_t cycle = 0ul;
    uint64_t reset_mask = UINT64_ONE_MASK;
    uint64_t init_idx = last_bit_idx & UINT64_BIT_MASK, cur_idx = init_idx;

    do {
        reset_mask &= get<1>(index_reset_mask(cur_idx));
        cur_idx = cur_idx + sieve_step;

        if (cur_idx >= UINT64_BITS)    {
            masks.push_back(reset_mask);
            reset_mask = UINT64_ONE_MASK;
            cur_idx &= UINT64_BIT_MASK;
        }

        if (cur_idx == init_idx)
            ++cycle;
    } while (cycle < repeats);

    assert(masks.size() == prime_number * repeats);

    return masks;
}

constexpr array<uint64_t, 8> compute_wheel_2_3_5()    {
    uint64_t wheel_idx = 0ul;
    array<uint64_t, 8ul> wheel = { };
    constexpr uint64_t wheel_size = 2ul * 3ul * 5ul;
    for (uint64_t n = 0ul; n < wheel_size; ++n)   {
        bool div_by_2 = (n % 2ul) == 0ul;
        bool div_by_3 = (n % 3ul) == 0ul;
        bool div_by_5 = (n % 5ul) == 0ul;
        bool div_by_2_3_5 = div_by_2 || div_by_3 || div_by_5;
        if (!div_by_2_3_5)
            wheel[wheel_idx++] = n;
    }

    return wheel;
}

constexpr array<int8_t, 30> compute_modulo_to_wheel_idx()   {
    constexpr uint64_t wheel_size = 2ul * 3ul * 5ul;
    array<int8_t, wheel_size> mod_to_idx = { };
    for (uint64_t n = 0; n < wheel_size; ++n)
        mod_to_idx[n] = -1;

    constexpr array<uint64_t, 8> wheel = compute_wheel_2_3_5();
    for (uint64_t r = 0; r < wheel.size(); ++r)
        mod_to_idx[wheel[r]] = r;

    return mod_to_idx;
}

const array<uint64_t, 8> Eratosthenes::wheel = compute_wheel_2_3_5();
const array<int8_t, 30> Eratosthenes::modulo_to_idx = compute_modulo_to_wheel_idx();
