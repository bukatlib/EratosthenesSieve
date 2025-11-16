#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include "cpu_info.h"
#include "eratosthenes_sieve.h"

#ifdef WHEEL_2_3_5
#define STEPS(prime_steps)  \
        prime_steps.step0, prime_steps.step1, prime_steps.step2, prime_steps.step3, \
        prime_steps.step4, prime_steps.step5, prime_steps.step6, prime_steps.step7
#define STEP_IN_BITS(steps) \
        steps[0] + steps[1] + steps[2] + steps[3] + steps[4] + steps[5] + steps[6] + steps[7]
#endif

#ifdef WHEEL_2_3
#define STEPS(prime_steps) prime_steps.step0, prime_steps.step1
#define STEP_IN_BITS(steps) steps[0] + steps[1]
#endif

#ifdef WHEEL_2
#define STEPS(prime_steps) prime_steps.step0
#define STEP_IN_BITS(steps) steps[0]
#endif


using namespace std;
using namespace std::chrono;

Eratosthenes::Eratosthenes(uint64_t primes_to) : 
    sieve_bits(sieve_size_bits(primes_to)), sieve_size(sieve_size_elems(sieve_bits)),
    primes_to(primes_to), sieve(new uint64_t[sieve_size]), segment_size(DEFAULT_SEGMENT_SIZE)
{
    #pragma omp parallel for if(sieve_bits > PARALLEL_SIEVE_MIN_BITS) num_threads(INIT_THREADS)
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
    // Segment size should be a multiple of uint64_t size!
    segment_size *= BITS_PER_BYTE;
    assert(segment_size % UINT64_BITS == 0);

    // Setup the minimal segment size to avoid reseting the previous bits.
    segment_size = max(segment_size, WHEEL_MIN_SEGMENT_SIZE);
}

Eratosthenes::~Eratosthenes()	{
	delete[] sieve;
}

uint64_t Eratosthenes::bit_size() const {
    return sieve_bits;
}

uint64_t Eratosthenes::bytes_allocated() const {
    return UINT64_BYTES * sieve_size;
}

void Eratosthenes::sieve_primes()    {
    // Handle the case when sieve is empty.
    if (sieve_bits == 0ul)
        return;

    uint64_t prime_bit_idx = 0;
    uint64_t segment_start = 0ul, segment_end = min(segment_size, sieve_bits);

    prime_bit_masks bit_masks;
    prime_wheel_data wheel_data;

    while ((prime_bit_idx = find_next_set_bit(prime_bit_idx)) < sieve_bits)  {
        if (prime_bit_idx >= segment_end)   {
            segment_start = segment_end;
            segment_end = min(segment_start + segment_size, sieve_bits);
            sieve_bit_range(bit_masks, wheel_data, segment_start, segment_end);
        }

        assert(prime_bit_idx < segment_end);

        // The newly found prime in the sieve image.
        // Seed only with the primes up to the root square of n.
        uint64_t prime = sieve_bit_to_number(prime_bit_idx);
        if (prime * prime > primes_to)  {
            break;
        }

        // Sieve multiples of the prime, primes larger or equal to prime^2 are considered.
        prime_wheel_steps prime_steps = wheel_steps(prime);
        sieve_bit_range_medium(prime_steps, segment_end);

        if (prime <= min(MAX_SMALL_PRIME, UINT64_BITS))
            bit_masks.push_back(seeding_pattern(prime));
        else
            wheel_data.push_back(prime_steps);

        prime_bit_idx++;
    }

    uint64_t chunk_size = WORKCHUNK_SEGMENTS * segment_size;
    uint64_t chunk_last_idx = segment_end + chunk_size * ((sieve_bits - segment_end) / chunk_size);

    #pragma omp parallel if(sieve_bits >= PARALLEL_SIEVE_MIN_BITS)
    {
        #pragma omp for schedule(dynamic)
        for (uint64_t b = segment_end; b < chunk_last_idx; b += chunk_size)   {
            uint64_t start_bit = b, end_bit = b + chunk_size;
            sieve_bit_range(bit_masks, wheel_data, start_bit, end_bit);
        }

        #pragma omp for schedule(static)
        for (uint64_t b = chunk_last_idx; b < sieve_bits; b += segment_size)    {
            uint64_t segment_start = b, segment_end = min(sieve_bits, b + segment_size);
            sieve_bit_range(bit_masks, wheel_data, segment_start, segment_end);
        }
    }
}

uint64_t Eratosthenes::prime_count() const {
    if (primes_to >= WHEEL_MIN_PRIME)  {
        uint64_t size_aligned = sieve_size & ALIGN_FOUR_MASK;
        uint64_t count0 = 0ul, count1 = 0ul, count2 = 0ul, count3 = 0ul, count4 = 0ul;

        #pragma omp parallel for reduction(+: count0,count1,count2,count3) if(sieve_bits > PARALLEL_SIEVE_MIN_BITS) num_threads(INIT_THREADS)
        for (uint64_t uint64_idx = 0ul; uint64_idx < size_aligned; uint64_idx += 4ul)  {
            count0 += popcount(sieve[uint64_idx + 0ul]);
            count1 += popcount(sieve[uint64_idx + 1ul]);
            count2 += popcount(sieve[uint64_idx + 2ul]);
            count3 += popcount(sieve[uint64_idx + 3ul]);
        }

        for (uint64_t uint64_idx = size_aligned; uint64_idx < sieve_size; ++uint64_idx)
            count4 += popcount(sieve[uint64_idx]);

        // Add implicit primes based on the wheel type, for example, primes 2, 3, 5 for 2x3x5 wheel.
        return WHEEL_IMPLICIT_PRIMES_COUNT + count0 + count1 + count2 + count3 + count4;
    } else {
        // primes to:                                         0    1    2    3    4    5    6
        static constexpr array<uint64_t, 7ul> prime2count = { 0ul, 0ul, 1ul, 2ul, 2ul, 3ul, 3ul };
        assert(primes_to < prime2count.size());
        return prime2count[primes_to];
    }
}

bool Eratosthenes::is_prime(uint64_t prime_candidate) const {
    static constexpr array<uint64_t, WHEEL_IMPLICIT_PRIMES_COUNT> wheel_implicit_primes = WHEEL_IMPLICIT_PRIMES;
    if (count(wheel_implicit_primes.cbegin(), wheel_implicit_primes.cend(), prime_candidate) == 1ul)    {
        // It is a prime used in the wheel's basis.
        return true;
    } else if (prime_candidate <= primes_to)    {
        uint64_t bit_idx = 0ul;
        if (is_in_image(prime_candidate, &bit_idx)) {
            // Number stored in the sieve binary image, check the bit.
            auto [idx, mask] = index_mask(bit_idx);
            return (sieve[idx] & mask) > 0ul;
        } else {
            // Number not in the sieve image, a multiple of implicit prime(s).
            return false;
        }
    } else {
        throw out_of_range("Increase the size of the sieve, too large prime to check!");
    }
}

vector<uint64_t> Eratosthenes::primes() const {
    vector<uint64_t> primes;
    primes.reserve(prime_count());

    for (uint64_t implicit_prime: WHEEL_IMPLICIT_PRIMES) {
        if (implicit_prime <= primes_to)
            primes.push_back(implicit_prime);
    }

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

    for (uint64_t implicit_prime: WHEEL_IMPLICIT_PRIMES) {
        if (implicit_prime <= primes_to)
            OUT<<implicit_prime<<"\n";
    }

    if (sieve_size > 0ul)  {
        array<uint64_t, 4ul> prime_buffer;
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

tuple<uint64_t, uint64_t> Eratosthenes::index_mask(uint64_t bit_idx)   {
    // Index to sieve vector and set/get bitmask.
    return {bit_idx >> UINT64_IDX_SHIFT, 1ul<<(UINT64_BIT_MASK & bit_idx)};
}

tuple<uint64_t, uint64_t> Eratosthenes::index_reset_mask(uint64_t bit_idx)   {
    // Index to sieve vector and set/get bitmask.
    return {bit_idx >> UINT64_IDX_SHIFT, rotl(CLEAR_BIT0_MASK, UINT64_BIT_MASK & bit_idx)};
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

bool Eratosthenes::is_in_image(uint64_t number, uint64_t* bit_idx) {
    if ((number & 0x1) == 0ul)  {
        // Divisible by 2 -> not in the sieve image for wheels 2, 2x3, and 2x3x5.
        return false;
    } else {
        uint64_t wheel_count = number / WHEEL_CIRCUMFERENCE;
        uint64_t modulo = number - WHEEL_CIRCUMFERENCE * wheel_count;
        int8_t wheel_bit_idx = modulo_to_idx[modulo];
        if (wheel_bit_idx >= 0)   {
            // Not divisible by any implicit prime -> in the sieve binary image.
            if (bit_idx != nullptr)
                *bit_idx = WHEEL_STEPS * wheel_count + wheel_bit_idx;
            return true;
        } else {
            // Divisible by an implicit prime -> not in the sieve image.
            return false;
        }
    }
}

uint64_t Eratosthenes::sieve_bit_to_number(uint64_t bit_idx)   {
    // Calculate the wheel index in the sieve and the modulo.
    uint64_t wheel_count = bit_idx / WHEEL_STEPS;
    uint64_t wheel_bit_idx = bit_idx - WHEEL_STEPS * wheel_count;
    uint64_t number = WHEEL_CIRCUMFERENCE * wheel_count + wheel[wheel_bit_idx];
    return number;
}

uint64_t Eratosthenes::sieve_size_bits(uint64_t primes_to)  {
    // Count the number of bits needed for the whole wheels first.
    uint64_t wheel_count = primes_to / WHEEL_CIRCUMFERENCE;
    uint64_t bitcount_whole_wheels = WHEEL_STEPS * wheel_count;

    // Define the modulo to determine the number of required bits for the partial (last) wheel.
    uint64_t bitcount_total = bitcount_whole_wheels, idx = 0ul;
    uint64_t modulo_circumference = primes_to - WHEEL_CIRCUMFERENCE * wheel_count;
    while (idx < wheel.size() && modulo_circumference >= wheel[idx++])
        bitcount_total++;

    return bitcount_total;
}

uint64_t Eratosthenes::sieve_size_elems(uint64_t sieve_bits)  {
    // Round up the size to 64 bit multiple.
    return (sieve_bits + UINT64_BITS - 1) / UINT64_BITS;
}

#ifdef WHEEL_2_3_5
constexpr array<uint64_t, WHEEL_STEPS> Eratosthenes::compute_wheel_2_3_5()    {
    uint64_t wheel_idx = 0ul;
    array<uint64_t, WHEEL_STEPS> wheel = { };
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

constexpr array<int8_t, WHEEL_CIRCUMFERENCE> Eratosthenes::compute_modulo_to_wheel_idx()   {
    array<int8_t, WHEEL_CIRCUMFERENCE> mod_to_idx = { };
    for (uint64_t n = 0; n < WHEEL_CIRCUMFERENCE; ++n)
        mod_to_idx[n] = -1;

    constexpr auto wheel = compute_wheel_2_3_5();
    for (uint64_t r = 0; r < wheel.size(); ++r)
        mod_to_idx[wheel[r]] = r;

    return mod_to_idx;
}

alignas(64) const array<uint64_t, WHEEL_STEPS> Eratosthenes::wheel = Eratosthenes::compute_wheel_2_3_5();
alignas(32) const array<int8_t, WHEEL_CIRCUMFERENCE> Eratosthenes::modulo_to_idx = Eratosthenes::compute_modulo_to_wheel_idx();
#endif

#ifdef WHEEL_2_3
alignas(64) const array<uint64_t, WHEEL_STEPS> Eratosthenes::wheel = { 1ul, 5ul };
alignas(32) const array<int8_t, WHEEL_CIRCUMFERENCE> Eratosthenes::modulo_to_idx = { -1, 0, -1, -1, -1, 1 };
#endif

#ifdef WHEEL_2
alignas(64) const array<uint64_t, WHEEL_STEPS> Eratosthenes::wheel = { 1ul };
alignas(32) const array<int8_t, WHEEL_CIRCUMFERENCE> Eratosthenes::modulo_to_idx = { -1, 0 };
#endif

Eratosthenes::prime_wheel_steps Eratosthenes::wheel_steps(uint64_t prime)   {
    // Start from the square of the found prime.
    uint64_t current_number = prime * prime;

    array<bool, WHEEL_STEPS> modulo_idx_found = {};
    array<uint64_t, WHEEL_STEPS> wheel_modulo_idxs = {};

    do {
        uint64_t bit_idx = 0ul;
        if (is_in_image(current_number, &bit_idx))    {
            #ifdef WHEEL_2_3_5
            uint64_t modulo_idx = bit_idx & 0x7;
            #else
            uint64_t modulo_idx = bit_idx - WHEEL_STEPS * (current_number / WHEEL_CIRCUMFERENCE);
            #endif
            assert(modulo_idx_found[modulo_idx] == false);
            wheel_modulo_idxs[modulo_idx] = bit_idx;
            modulo_idx_found[modulo_idx] = true;
        }

        current_number += prime;
    } while (find(modulo_idx_found.cbegin(), modulo_idx_found.cend(), false) != modulo_idx_found.cend());

    array<uint64_t, WHEEL_STEPS + 1ul> bit_idxs_shifted = {};
    auto min_bit_idx_sit = min_element(wheel_modulo_idxs.cbegin(), wheel_modulo_idxs.cend());
    for (uint64_t i = 0ul; i < wheel_modulo_idxs.size(); ++i)
        bit_idxs_shifted[i] = wheel_modulo_idxs[i] - *min_bit_idx_sit;
    bit_idxs_shifted[WHEEL_STEPS] = WHEEL_STEPS * prime;

    array<uint64_t, WHEEL_STEPS> prime_steps = {};
    sort(bit_idxs_shifted.begin(), bit_idxs_shifted.end());
    for (uint64_t i = 0ul; i + 1ul < bit_idxs_shifted.size(); ++i)
        prime_steps[i] = bit_idxs_shifted[i + 1ul] - bit_idxs_shifted[i];

    assert(accumulate(prime_steps.cbegin(), prime_steps.cend(), 0ul) == WHEEL_STEPS * prime);

    prime_wheel_steps wheel_steps_struct;

    #ifdef WHEEL_2_3_5
    wheel_steps_struct.step_idx = 0;
    wheel_steps_struct.bit_idx = *min_bit_idx_sit;
    wheel_steps_struct.step0 = prime_steps[0];
    wheel_steps_struct.step1 = prime_steps[1];
    wheel_steps_struct.step2 = prime_steps[2];
    wheel_steps_struct.step3 = prime_steps[3];
    wheel_steps_struct.step4 = prime_steps[4];
    wheel_steps_struct.step5 = prime_steps[5];
    wheel_steps_struct.step6 = prime_steps[6];
    wheel_steps_struct.step7 = prime_steps[7];
    #endif

    #ifdef WHEEL_2_3
    wheel_steps_struct.step_idx = 0;
    wheel_steps_struct.bit_idx = *min_bit_idx_sit;
    wheel_steps_struct.step0 = prime_steps[0];
    wheel_steps_struct.step1 = prime_steps[1];
    #endif

    #ifdef WHEEL_2
    wheel_steps_struct.bit_idx = *min_bit_idx_sit;
    wheel_steps_struct.step0 = prime_steps[0];
    #endif

    return wheel_steps_struct;
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
                // We allow the holes up to one uint64_t; no bigger holes are allowed.
                assert(uint64_idx - pattern_idx <= 2);

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

    for (uint64_t n = 0; n < repeats; ++n)  {
        for (uint64_t mask: cycle_masks)
            masks.push_back(mask);
    }

    return masks;
}

void Eratosthenes::sieve_bit_range(const prime_bit_masks& masks, const prime_wheel_data& wheel_data,
                                    uint64_t start_bit, uint64_t end_bit)    {

    vector<uint64_t> pattern_idxs(masks.size(), 0ul);

    for (uint64_t i = 0ul; i < masks.size(); ++i)  {
        assert(start_bit >= segment_size);
        uint64_t mask_count = masks[i].size();
        assert(mask_count > 0ul);
        assert((start_bit & UINT64_BIT_MASK) == 0ul);
        uint64_t pattern_idx = (start_bit / UINT64_BITS) % mask_count;
        pattern_idxs[i] = pattern_idx;
    }

    prime_wheel_data wheel_steps;
    wheel_steps.reserve(wheel_data.size());
    uint64_t min_idx_big_step = wheel_data.size();
    constexpr uint64_t MAX_STEPS_FOR_SMALL = 2ul;

    for (uint64_t i = 0ul; i < wheel_data.size(); ++i)  {
        prime_wheel_steps prime_data = wheel_data[i];
        uint32_t steps[WHEEL_STEPS] = { STEPS(prime_data) };
        uint64_t sieve_step_in_bits = STEP_IN_BITS(steps);

        uint64_t init_bit = prime_data.bit_idx;
        uint64_t step_mult = (max(start_bit, init_bit) - init_bit) / sieve_step_in_bits;
        uint64_t start_bit_wheel_idx = init_bit + step_mult * sieve_step_in_bits;
        assert((start_bit_wheel_idx - init_bit) % sieve_step_in_bits == 0ul);

        #ifdef WHEEL_2
        if (start_bit_wheel_idx < start_bit)
            start_bit_wheel_idx += sieve_step_in_bits;
        #else
        uint8_t step_idx = prime_data.step_idx;
        while (start_bit_wheel_idx < start_bit) {
            start_bit_wheel_idx += steps[step_idx];
            step_idx = (step_idx + 1u) & WHEEL_STEPS_MASK;
        }
        #endif

        if (start_bit_wheel_idx < end_bit)   {
            #if not defined WHEEL_2
            prime_data.step_idx = step_idx;
            #endif
            prime_data.bit_idx = start_bit_wheel_idx;
            if (MAX_STEPS_FOR_SMALL*sieve_step_in_bits > segment_size)
                min_idx_big_step = min(min_idx_big_step, i);
            wheel_steps.push_back(prime_data);
        }
    }

    min_idx_big_step = min(min_idx_big_step, wheel_steps.size());

    for (uint64_t segment_start = start_bit; segment_start < end_bit; segment_start += segment_size)    {
        uint64_t segment_end = min(end_bit, segment_start + segment_size); 

        // Zero multiples of small primes by using bit masks.
        sieve_bit_range_small(masks, pattern_idxs, segment_start, segment_end);

        // Reset idividual bits for medium primes.
        for (uint64_t j = 0ul; j < min_idx_big_step; ++j)
            sieve_bit_range_medium(wheel_steps[j], segment_end);

        // Reset idividual bits for large primes.
        for (uint64_t j = min_idx_big_step; j < wheel_steps.size(); ++j)
            sieve_bit_range_large(wheel_steps[j], segment_end);
    }
}

void Eratosthenes::sieve_bit_range_small(const prime_bit_masks& masks, vector<uint64_t>& pattern_idxs, uint64_t start_bit, uint64_t end_bit)   {
    // Use pattern to seed small primes.
    for (uint64_t i = 0ul; i < masks.size(); ++i) {
        const vector<uint64_t>& pattern = masks[i];
        uint64_t pattern_chunks = pattern.size();
        uint64_t segment_bits = end_bit - start_bit;
        uint64_t partial_bits = segment_bits & UINT64_BIT_MASK;
        uint64_t sieve_idx = (start_bit>>UINT64_IDX_SHIFT), pattern_idx = pattern_idxs[i];
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
}

#ifdef WHEEL_2_3_5
void Eratosthenes::sieve_bit_range_medium(prime_wheel_steps& prime_steps, uint64_t end_bit)    {
    uint64_t start_bit = prime_steps.bit_idx;
    if (start_bit >= end_bit)
        return;

    uint8_t step_idx = prime_steps.step_idx;
    uint32_t steps[WHEEL_STEPS] = { STEPS(prime_steps) };

    uint64_t offset0 = 0ul;
    uint64_t offset1 = offset0 + steps[(step_idx + 0u) & WHEEL_STEPS_MASK];
    uint64_t offset2 = offset1 + steps[(step_idx + 1u) & WHEEL_STEPS_MASK];
    uint64_t offset3 = offset2 + steps[(step_idx + 2u) & WHEEL_STEPS_MASK];
    uint64_t offset4 = offset3 + steps[(step_idx + 3u) & WHEEL_STEPS_MASK];
    uint64_t offset5 = offset4 + steps[(step_idx + 4u) & WHEEL_STEPS_MASK];
    uint64_t offset6 = offset5 + steps[(step_idx + 5u) & WHEEL_STEPS_MASK];
    uint64_t offset7 = offset6 + steps[(step_idx + 6u) & WHEEL_STEPS_MASK];
    uint64_t prime_step_in_bits = offset7 + steps[(step_idx + 7u) & WHEEL_STEPS_MASK];
    uint64_t max_unroll_idx = end_bit - min(end_bit, prime_step_in_bits);

    while (start_bit < max_unroll_idx)  {
        reset_bit(start_bit + offset0);
        reset_bit(start_bit + offset1);
        reset_bit(start_bit + offset2);
        reset_bit(start_bit + offset3);
        reset_bit(start_bit + offset4);
        reset_bit(start_bit + offset5);
        reset_bit(start_bit + offset6);
        reset_bit(start_bit + offset7);
        start_bit += prime_step_in_bits;
    }

    assert(start_bit + prime_step_in_bits >= end_bit);

    uint8_t next_step_idx = step_idx;
    while (start_bit < end_bit)  {
        reset_bit(start_bit);
        start_bit += steps[next_step_idx];
        next_step_idx = (next_step_idx + 1u) & WHEEL_STEPS_MASK;
    }

    prime_steps.step_idx = next_step_idx;
    prime_steps.bit_idx = start_bit;
}
#endif

#ifdef WHEEL_2_3
void Eratosthenes::sieve_bit_range_medium(prime_wheel_steps& prime_steps, uint64_t end_bit)    {
    uint64_t start_bit = prime_steps.bit_idx;
    if (start_bit >= end_bit)
        return;

    uint8_t step_idx = prime_steps.step_idx;
    uint32_t steps[WHEEL_STEPS] = { STEPS(prime_steps) };

    uint64_t offset0 = 0ul;
    uint64_t offset1 = offset0 + steps[(step_idx + 0u) & WHEEL_STEPS_MASK];
    uint64_t prime_step_in_bits = offset1 + steps[(step_idx + 1u) & WHEEL_STEPS_MASK];
    uint64_t offset2 = offset0 + prime_step_in_bits;
    uint64_t offset3 = offset1 + prime_step_in_bits;
    uint64_t loop_step = 2ul * prime_step_in_bits;
    uint64_t max_unroll_idx = end_bit - min(end_bit, loop_step);

    while (start_bit < max_unroll_idx)  {
        reset_bit(start_bit + offset0);
        reset_bit(start_bit + offset1);
        reset_bit(start_bit + offset2);
        reset_bit(start_bit + offset3);
        start_bit += loop_step;
    }

    assert(start_bit + loop_step >= end_bit);

    uint8_t next_step_idx = step_idx;
    while (start_bit < end_bit)  {
        reset_bit(start_bit);
        start_bit += steps[next_step_idx];
        next_step_idx = (next_step_idx + 1u) & WHEEL_STEPS_MASK;
    }

    prime_steps.step_idx = next_step_idx;
    prime_steps.bit_idx = start_bit;
}
#endif

#ifdef WHEEL_2
void Eratosthenes::sieve_bit_range_medium(prime_wheel_steps& prime_steps, uint64_t end_bit)    {
    uint64_t start_bit = prime_steps.bit_idx;
    if (start_bit >= end_bit)
        return;

    uint64_t offset0 = 0ul, offset1 = prime_steps.step0;
    uint64_t offset2 = 2ul * offset1, offset3 = offset2 + offset1;
    uint64_t loop_step = 4ul * offset1;
    uint64_t max_unroll_idx = end_bit - min(end_bit, loop_step);

    while (start_bit < max_unroll_idx)  {
        reset_bit(start_bit + offset0);
        reset_bit(start_bit + offset1);
        reset_bit(start_bit + offset2);
        reset_bit(start_bit + offset3);
        start_bit += loop_step;
    }

    assert(start_bit + loop_step >= end_bit);

    while (start_bit < end_bit)  {
        reset_bit(start_bit);
        start_bit += offset1;
    }

    prime_steps.bit_idx = start_bit;
}
#endif

#if defined WHEEL_2_3_5 || defined WHEEL_2_3
void Eratosthenes::sieve_bit_range_large(prime_wheel_steps& prime_steps, uint64_t end_bit)    {
    uint64_t start_bit = prime_steps.bit_idx;
    if (start_bit >= end_bit)
        return;

    uint32_t steps[WHEEL_STEPS] = { STEPS(prime_steps) };

    uint8_t next_step_idx = prime_steps.step_idx;
    while (start_bit < end_bit)  {
        reset_bit(start_bit);
        start_bit += steps[next_step_idx];
        next_step_idx = (next_step_idx + 1) & WHEEL_STEPS_MASK;
    }

    prime_steps.step_idx = next_step_idx;
    prime_steps.bit_idx = start_bit;
}
#endif

#ifdef WHEEL_2
void Eratosthenes::sieve_bit_range_large(prime_wheel_steps& prime_steps, uint64_t end_bit)    {
    uint64_t start_bit = prime_steps.bit_idx;
    if (start_bit >= end_bit)
        return;

    uint64_t step = prime_steps.step0;
    while (start_bit < end_bit)  {
        reset_bit(start_bit);
        start_bit += step;
    }

    prime_steps.bit_idx = start_bit;
}
#endif
