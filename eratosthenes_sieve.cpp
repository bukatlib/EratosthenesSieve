#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "cpu_info.h"
#include "eratosthenes_sieve.h"

// TODO: remove
#include <iomanip>

using namespace std;
using namespace std::chrono;

#pragma pack(1)
struct test_bitfield {
    uint64_t step_idx : 3;
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

    // Clear number '1' in the sieve image as it is not a prime.
    if (sieve_size > 0ul)   {
        constexpr uint64_t number_1_clear_mask = UINT64_ONE_MASK - 1ul;
        sieve[0] &= number_1_clear_mask;
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
    // TODO: update
    //assert(segment_size >= 2048ul);
    //segment_size = 1024ul;
    // segment_size = 384*1024ul;
    // Segment size should be a multiple of uint64_t size!
    assert(segment_size % UINT64_BITS == 0);
}

uint64_t test(uint64_t sieve_bits) {
    test_bitfield tf;
    tf.step_idx = (sieve_bits & 0x1ul) + 1u;
    tf.bit_idx = sieve_bits;
    tf.step0 = sieve_bits>>0; 
    tf.step1 = sieve_bits>>1; 
    tf.step2 = sieve_bits>>2; 
    tf.step3 = sieve_bits>>3; 
    tf.step4 = sieve_bits>>4; 
    tf.step5 = sieve_bits>>5; 
    tf.step6 = sieve_bits>>6; 
    tf.step7 = sieve_bits>>7; 
    return tf.step_idx * tf.bit_idx + tf.step0 + tf.step1 + tf.step2 + tf.step3 + tf.step4 + tf.step5 + tf.step6 + tf.step7;
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

        // The newly found prime in the sieve image.
        // Seed only with the primes up to the root square of n.
        uint64_t prime = sieve_bit_to_number(prime_bit_idx);
        if (prime * prime > primes_to)  {
            break;
        }

        // Sieve multiples of the prime, primes larger or equal to prime^2 are considered.
        array<uint64_t, 8> wheel_modulo_idxs = Eratosthenes::wheel_modulo_init_idxs(prime);
        for (uint64_t start_bit_idx: wheel_modulo_idxs) {
            uint64_t bit_idx = start_bit_idx;
            while (bit_idx < segment_end)   {
                reset_bit(bit_idx);
                bit_idx += 8ul * prime;
            }
        }

        if (prime <= min(MAX_SMALL_PRIME, UINT64_BITS))   {
            primes_meta.push_back({prime, wheel_modulo_idxs, seeding_pattern(prime)});
            ++patterns;
        } else {
            primes_meta.push_back({prime, wheel_modulo_idxs, {}});
        }
        
        prime_bit_idx++;
    }

    /*
    uint64_t total_size = 0ul;
    for (const prime_seed_metainfo& meta: primes_meta)
        total_size += 8ul * meta.reset_masks.size();
    cout<<"Total patterns: "<<patterns<<" ("<<total_size<<" bytes)"<<endl;
    */
    uint64_t val = test(sieve_bits);
    cout<<"val: "<<val<<endl;
    return;

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
    if (primes_to >= 7ul)  {
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

        // Add three for primes 2, 3, and 5.
        found_primes = 3ul + count0 + count1 + count2 + count3 + count4;
    } else if (primes_to >= 5ul) {
        // Primes 2, 3, 5 but not 7 and greater.
        found_primes = 3ul;
    } else if (primes_to >= 3ul) {
        // Primes 2, 3 but not 5 and greater.
        found_primes = 2ul;
    } else if (primes_to >= 2ul) {
        // Prime 2 but not 3 and greater.
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
    if (primes_to >= 3ul)
        primes.push_back(3ul);
    if (primes_to >= 5ul)
        primes.push_back(5ul);
    
    if (sieve_size > 0ul)  {
        uint64_t prime_bit_idx = 0;
        while ((prime_bit_idx = find_next_set_bit(prime_bit_idx)) < sieve_bits) {
            primes.push_back(sieve_bit_to_number(prime_bit_idx));
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
    if (primes_to >= 3ul)
        OUT<<3<<"\n";
    if (primes_to >= 5ul)
        OUT<<5<<"\n";

    if (sieve_size > 0ul)  {
        array<uint64_t, 4> prime_buffer;
        uint64_t buffer_idx = 0ul, prime_bit_idx = 0ul;
        while ((prime_bit_idx = find_next_set_bit(prime_bit_idx)) < sieve_bits) {
            uint64_t prime_val = sieve_bit_to_number(prime_bit_idx);

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

uint64_t Eratosthenes::sieve_bit_to_number(uint64_t bit_idx)   {
    // Calculate the wheel 2x3x5 index in the sieve and the modulo.
    constexpr uint64_t basis_size = 8ul;
    constexpr uint64_t wheel_circumference = 30ul;
    uint64_t wheel_count = bit_idx / basis_size;
    uint64_t wheel_bit_idx = bit_idx - basis_size * wheel_count;
    uint64_t number = wheel_circumference * wheel_count + wheel[wheel_bit_idx];
    return number;
}

bool Eratosthenes::is_in_image(uint64_t number, uint64_t* bit_idx) {
    if ((number & 0x1) == 0ul)  {
        // Divisible by 2 -> not in the sieve image.
        return false;
    } else {
        constexpr uint64_t wheel_circumference = 30ul;
        uint64_t wheel_count = number / wheel_circumference;
        uint64_t modulo = number - wheel_circumference * wheel_count;
        int8_t wheel_bit_idx = Eratosthenes::modulo_to_idx[modulo];
        if (wheel_bit_idx >= 0)   {
            // Not divisible by 2, 3, or 5 -> in the sieve binary image.
            if (bit_idx != nullptr)
                *bit_idx = 8ul * wheel_count + wheel_bit_idx;
            return true;
        } else {
            // Divisible by 3 or 5 -> not in the sieve image.
            return false;
        }
    }
}

array<uint64_t, 8> Eratosthenes::wheel_modulo_init_idxs(uint64_t prime)   {
    // Start from the square of the found prime.
    uint64_t current_number = prime * prime;

    array<bool, 8> modulo_idx_found = { };
    array<uint64_t, 8> wheel_modulo_idxs = { };

    do {
        uint64_t bit_idx = 0ul;
        if (Eratosthenes::is_in_image(current_number, &bit_idx))    {
            uint64_t byte_bit_idx = bit_idx & 0x7;
            assert(modulo_idx_found[byte_bit_idx] == false);
            wheel_modulo_idxs[byte_bit_idx] = bit_idx;
            modulo_idx_found[byte_bit_idx] = true;
        }

        current_number += prime;
    } while (find(modulo_idx_found.cbegin(), modulo_idx_found.cend(), false) != modulo_idx_found.cend());

    array<uint64_t, 8> bit_steps = { };
    uint64_t min_bit = *min_element(wheel_modulo_idxs.cbegin(), wheel_modulo_idxs.cend());
    for (uint64_t i = 0; i < 8ul; ++i)
        bit_steps[i] = wheel_modulo_idxs[i] - min_bit;

    cout<<"prime "<<prime<<": "<<*max_element(bit_steps.cbegin(), bit_steps.cend())<<endl;

    return wheel_modulo_idxs;
}

uint64_t Eratosthenes::sieve_size_bits(uint64_t primes_to)  {
    // Count the number of bits needed for the whole wheels first.
    constexpr uint64_t wheel_bits = 8ul;
    constexpr uint64_t wheel_circumference = 30ul;
    uint64_t wheel_count = primes_to / wheel_circumference;
    uint64_t bitcount_whole_wheels = wheel_bits * wheel_count;

    // Define the modulo to determine the number of required bits for the partial (last) wheel.
    uint64_t bitcount_total = bitcount_whole_wheels, idx = 0ul;
    uint64_t modulo_circumference = primes_to - wheel_circumference * wheel_count;
    while (idx < Eratosthenes::wheel.size() && modulo_circumference >= Eratosthenes::wheel[idx++])
        bitcount_total++;

    return bitcount_total;
}

uint64_t Eratosthenes::sieve_size_bytes_alloc(uint64_t sieve_bits)  {
    // Round up the size to 64 bit multiple.
    return UINT64_BYTES * ((sieve_bits + UINT64_BITS - 1) / UINT64_BITS);
}

void Eratosthenes::seed_bit_range(const vector<prime_seed_metainfo>& primes_meta,
                    uint64_t num_of_patterns, uint64_t start_bit, uint64_t end_bit)    {

    vector<uint64_t> pattern_idxs(num_of_patterns, 0ul);

    for (uint64_t i = 0ul; i < num_of_patterns; ++i)  {
        assert(start_bit >= segment_size);
        uint64_t mask_count = primes_meta[i].reset_masks.size();
        assert(mask_count > 0ul);
        assert((start_bit & UINT64_BIT_MASK) == 0ul);
        uint64_t pattern_idx = (start_bit / UINT64_BITS) % mask_count;
        pattern_idxs[i] = pattern_idx;
    }

    vector<uint64_t> large_prime_idxs;
    vector<uint32_t> large_prime_steps;
    uint64_t prime_idx = num_of_patterns;
    assert(num_of_patterns <= primes_meta.size());
    uint64_t number_of_large_primes = primes_meta.size() - num_of_patterns;
    large_prime_idxs.reserve(8ul*number_of_large_primes);

    while (prime_idx < primes_meta.size())    {
        uint64_t sieve_step_in_bits = 8ul * primes_meta[prime_idx].prime;
        const array<uint64_t, 8ul>& modulo_init_idxs = primes_meta[prime_idx].modulo_init_idxs;

        for (uint64_t i = 0ul; i < modulo_init_idxs.size(); ++i)    {
            uint64_t init_modulo_bit = modulo_init_idxs[i];
            uint64_t start_bit_modulo = max(start_bit, init_modulo_bit);
            uint64_t step_mult = (start_bit_modulo - init_modulo_bit + sieve_step_in_bits - 1ul) / sieve_step_in_bits;
            uint64_t cur_modulo_idx = init_modulo_bit + step_mult * sieve_step_in_bits;
            assert((cur_modulo_idx - init_modulo_bit) % sieve_step_in_bits == 0ul);
            assert((sieve_bit_to_number(cur_modulo_idx) % primes_meta[prime_idx].prime) == 0ul);
            if (cur_modulo_idx < end_bit)   {
                large_prime_idxs.push_back(cur_modulo_idx);
                large_prime_steps.push_back(static_cast<uint32_t>(sieve_step_in_bits));
            }
        }

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
            uint64_t sieve_idx = (segment_start >> UINT64_IDX_SHIFT), pattern_idx = pattern_idxs[i];
            uint64_t pattern_writes = (segment_bits >> UINT64_IDX_SHIFT) + (partial_bits == 0ul ? 0ul : 1ul);
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

            if (pattern_idx + partial_pattern_count >= pattern_chunks)  {
                pattern_idxs[i] = pattern_writes;
            } else {
                pattern_idxs[i] = pattern_idx + partial_pattern_count;
            }
        }

        // Reset idividual bits for large primes.
        for (uint64_t i = 0; i < large_prime_idxs.size(); ++i)   {
            uint64_t bit_idx = large_prime_idxs[i];
            uint64_t one_step = large_prime_steps[i];

            while (bit_idx < segment_end)  {
                reset_bit(bit_idx);
                bit_idx += one_step;
            }

            large_prime_idxs[i] = bit_idx;
        }
    }
}

vector<uint64_t> Eratosthenes::seeding_pattern(uint64_t prime) {
    // Calculate seeding patterns for small prime numbers.
    vector<uint64_t> cycle_masks;
    cycle_masks.reserve(prime);

    uint64_t mask = UINT64_ONE_MASK;
    uint64_t pattern_idx = 0ul, current_number = prime;

    while (cycle_masks.size() < prime)    {
        uint64_t bit_idx = 0ul;
        if (is_in_image(current_number, &bit_idx))   {
            auto [uint64_idx, reset_bit_mask] = index_reset_mask(bit_idx);
            assert((sieve_bit_to_number(bit_idx) % prime) == 0ul);

            if (uint64_idx > pattern_idx)    {
                // We allow the holes up to one uint64_t.
                assert(uint64_idx - pattern_idx <= 2);  // No bigger holes in pattern allowed.

                cycle_masks.push_back(mask);
                if (pattern_idx + 1 < uint64_idx)   {
                    cycle_masks.push_back(UINT64_ONE_MASK);
                    ++pattern_idx;
                }

                mask = UINT64_ONE_MASK;
                ++pattern_idx;
            }

            mask &= reset_bit_mask;
        }

        current_number += prime;
    }

    vector<uint64_t> masks;
    uint64_t repeats = MAX_MASKS_PER_PRIME / prime;
    assert(repeats >= 1ul);
    masks.reserve(prime * repeats);

    //cout<<"prime "<<prime<<" with "<<repeats<<" (total "<<repeats * prime<<")"<<endl;

    for (uint64_t n = 0; n < repeats; ++n)  {
        for (uint64_t mask: cycle_masks)
            masks.push_back(mask);
    }

    /*
    cout<<"masks:";
    for (uint64_t idx = 0; idx < masks.size(); ++idx)   {
        if ((idx % prime) == 0ul)
            cout<<endl;
        cout<<" "<<hex<<masks[idx];
    }
    cout<<dec<<endl;
    */

    return masks;
}

constexpr array<uint64_t, Eratosthenes::WHEEL_NUMBER_OF_COPRIMES> Eratosthenes::compute_wheel_2_3_5()    {
    uint64_t wheel_idx = 0ul;
    array<uint64_t, WHEEL_NUMBER_OF_COPRIMES> wheel = { };
    for (uint64_t n = 0ul; n < WHEEL_CIRCUMFERENCE; ++n)   {
        bool div_by_2 = (n % 2ul) == 0ul;
        bool div_by_3 = (n % 3ul) == 0ul;
        bool div_by_5 = (n % 5ul) == 0ul;
        bool div_by_2_3_5 = div_by_2 || div_by_3 || div_by_5;
        if (!div_by_2_3_5)
            wheel[wheel_idx++] = n;
    }

    return wheel;
}

constexpr array<int8_t, Eratosthenes::WHEEL_CIRCUMFERENCE> Eratosthenes::compute_modulo_to_wheel_idx()   {
    array<int8_t, WHEEL_CIRCUMFERENCE> mod_to_idx = { };
    for (uint64_t n = 0; n < WHEEL_CIRCUMFERENCE; ++n)
        mod_to_idx[n] = -1;

    constexpr auto wheel = compute_wheel_2_3_5();
    for (uint64_t r = 0; r < wheel.size(); ++r)
        mod_to_idx[wheel[r]] = r;

    return mod_to_idx;
}

const array<uint64_t, Eratosthenes::WHEEL_NUMBER_OF_COPRIMES> Eratosthenes::wheel = Eratosthenes::compute_wheel_2_3_5();
const array<int8_t, Eratosthenes::WHEEL_CIRCUMFERENCE> Eratosthenes::modulo_to_idx = Eratosthenes::compute_modulo_to_wheel_idx();
