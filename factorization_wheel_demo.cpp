#include <algorithm>
#include <bit>
#include <map>
#include <iterator>
#include <iomanip>
#include <iostream>
#include <stdint.h>
#include <vector>
#include <string>

using namespace std;

struct candidate {
    uint64_t number;
    uint64_t modulo;
    uint64_t bit_idx;
};

bool div_by_2_3_5(uint64_t n)   {
    bool div_by_2 = (n % 2) == 0;
    bool div_by_3 = (n % 3) == 0;
    bool div_by_5 = (n % 5) == 0;
    return div_by_2 || div_by_3 || div_by_5;
}

int main() {
    const uint64_t factor = 2*3*5;

    map<uint64_t, uint64_t> modulo_to_bit = {
        {1, 0}, {7, 1}, {11, 2}, {13, 3}, {17, 4}, {19, 5}, {23, 6}, {29, 7}
    };

    map<uint64_t, uint64_t> bit_to_modulo;
    for (auto it = modulo_to_bit.cbegin(); it != modulo_to_bit.cend(); ++it)
        bit_to_modulo[it->second] = it->first;


    cout<<"BASE FOR 2x3x5 FACTORIZATION"<<endl;
    for (auto it = modulo_to_bit.cbegin(); it != modulo_to_bit.cend(); ++it)
        cout<<" "<<it->first;
    cout<<endl<<endl;

    vector<candidate> candidates;
    for (uint64_t n = 0; n < 10000; ++n)   {
        if (!div_by_2_3_5(n)) {
            uint64_t modulo = n % factor;
            candidates.emplace_back(n, modulo, modulo_to_bit[modulo]);
        }
    }

    cout<<"COPRIMES FOR 2x3x5 UP TO 1020"<<endl;
    for (uint64_t bit_idx = 0; bit_idx < 8; ++bit_idx)  {
        vector<candidate> modulo_candidates;
        uint64_t modulo = bit_to_modulo[bit_idx];

        auto filter = [=](const candidate& c) { return (c.number <= 1020) && (c.bit_idx == bit_idx); };
        copy_if(candidates.cbegin(), candidates.cend(), back_inserter(modulo_candidates), filter);

        cout<<setw(3)<<modulo<<":";
        for (const candidate& modulo_candidate: modulo_candidates)
            cout<<setw(5)<<modulo_candidate.number;
        cout<<endl;
    }

    for (uint64_t p: {7, 11, 13})   {
        cout<<endl<<"MULTIPLES OF "<<p<<" PRIME"<<endl;
        for (uint64_t bit_idx = 0; bit_idx < 8; ++bit_idx)  {
            vector<candidate> prime_mult_modulo;
            uint64_t modulo = bit_to_modulo[bit_idx];

            auto filter = [=](const candidate& c) { return (c.number <= 1020) && (c.bit_idx == bit_idx) && ((c.number % p) == 0);};
            copy_if(candidates.cbegin(), candidates.cend(), back_inserter(prime_mult_modulo), filter);

            cout<<setw(3)<<modulo<<":";
            for (const candidate& prime_mult: prime_mult_modulo)    {
                uint64_t byte_idx = prime_mult.number / factor;
                uint64_t bit_idx = 8ul * byte_idx + prime_mult.bit_idx;
                string info = to_string(byte_idx) + "/" + to_string(bit_idx);
                cout<<setw(8)<<info;
            }
            cout<<endl;
        }
    }

    for (uint64_t p: {7, 11, 13})   {
        map<uint64_t, uint64_t> sieve_masks;
        for (const candidate& c: candidates)    {
            if ((c.number % p) == 0)    {
                uint64_t byte_idx = c.number / factor;
                uint64_t bit_idx = 8ul * byte_idx + c.bit_idx;
                uint64_t arr_idx = byte_idx >> 3;
                uint64_t reset_mask = ~(1ul<<bit_idx);
                if (auto sit = sieve_masks.find(arr_idx); sit != sieve_masks.cend()) {
                    sit->second &= reset_mask;
                } else {
                    if (sieve_masks.size() < 3 * p)
                        sieve_masks[arr_idx] = reset_mask;
                    else
                        break;
                }
            }
        }

        uint64_t one_count = 0ul;
        cout<<endl<<"SIEVING MASKS FOR "<<p<<" PRIME"<<endl;
        for (auto it = sieve_masks.cbegin(); it != sieve_masks.cend(); ++it)  {
            cout<<hex<<" "<<it->second<<dec;
            one_count += popcount(it->second);
            if (((it->first + 1) % p) == 0) {
                cout<<endl;
            }
        }
        cout<<"... "<<(3 * p * 64 - one_count)<<" zeros of "<<3 * p * 64<<" total."<<endl;
    }

    return 0;
}
