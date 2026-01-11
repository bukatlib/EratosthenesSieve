#include <chrono>
#include <iostream>
#include <string>
#include <stdint.h>
#include "eratosthenes_sieve.h"

using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[])  {
	uint64_t to_value;
	string file_name;
    if (argc == 2)  {
        file_name = argv[1];
        if (!file_name.ends_with(".bin"))   {
            cerr<<argv[0]<<" <INPUT_IMAGE>"<<endl;
            return 1;
        }

        time_point<system_clock, duration<double>> tp1 = high_resolution_clock::now();

        Eratosthenes sieve(file_name);
        uint64_t bit_size = sieve.bit_size(), byte_size = sieve.bytes_allocated();
        double mbyte_size = byte_size / (1024.0 * 1024.0);
        cout<<"sieve size: "<<bit_size<<" bits ("<<byte_size<<" bytes, "<<mbyte_size<<" MBs)"<<endl;

        time_point<system_clock, duration<double>> tp2 = high_resolution_clock::now();
        double init_time = duration_cast<duration<double>>(tp2 - tp1).count();
        cout<<"Sieve up to "<<sieve.to_sieved()<<" loaded in "<<init_time<<" secs"<<endl;

        uint64_t prime_count = sieve.prime_count();

        time_point<system_clock, duration<double>> tp3 = high_resolution_clock::now();
        double count_time = duration_cast<duration<double>>(tp3 - tp2).count();
        double bandwidth = sieve.bytes_allocated() / (1000.0 * 1000.0 * count_time);
        cout<<"Found "<<prime_count<<" primes in "<<count_time<<" secs (bandwidth "<<bandwidth<<" MB/s)."<<endl;
    } else if (argc >= 3)  {
        string to_value_str = argv[1];
        string::size_type sz = 0; 
        to_value = stoull(to_value_str, &sz);
        string arg_output = argv[2];
        if (sz != to_value_str.size())  {
            cerr<<"Invalid '<PRIMES_TO>' argument!"<<endl;
            cerr<<argv[0]<<" <PRIMES_TO> <OUTPUT_FILE>"<<endl;
            return 1;
        }
        if (arg_output != "-" && !arg_output.empty())
            file_name = arg_output;

        time_point<system_clock, duration<double>> tp1 = high_resolution_clock::now();

        Eratosthenes sieve(to_value);
        uint64_t bit_size = sieve.bit_size(), byte_size = sieve.bytes_allocated();
        double mbyte_size = byte_size / (1024.0 * 1024.0);
        cout<<"sieve size: "<<bit_size<<" bits ("<<byte_size<<" bytes, "<<mbyte_size<<" MBs)"<<endl;

        time_point<system_clock, duration<double>> tp2 = high_resolution_clock::now();
        double init_time = duration_cast<duration<double>>(tp2 - tp1).count();
        cout<<"Sieve initialized in "<<init_time<<" secs"<<endl;

        sieve.sieve_primes();

        time_point<system_clock, duration<double>> tp3 = high_resolution_clock::now();
        double sieve_time = duration_cast<duration<double>>(tp3 - tp2).count();
        cout<<"Primes sieved in "<<sieve_time<<" secs."<<endl;

        uint64_t prime_count = sieve.prime_count();

        time_point<system_clock, duration<double>> tp4 = high_resolution_clock::now();
        double count_time = duration_cast<duration<double>>(tp4 - tp3).count();
        double bandwidth = sieve.bytes_allocated() / (1000.0 * 1000.0 * count_time);
        cout<<"Found "<<prime_count<<" primes in "<<count_time<<" secs (bandwidth "<<bandwidth<<" MB/s)."<<endl;

        if (!file_name.empty()) {
            if (file_name.ends_with(".bin"))
                sieve.write_image(file_name);
            else
                sieve.write_primes_to_file(file_name);

            time_point<system_clock, duration<double>> tp5 = high_resolution_clock::now();
            double write_time = duration_cast<duration<double>>(tp5 - tp4).count();
            cout<<"Primes written to '"<<file_name<<"' file in "<<write_time<<" secs."<<endl;
        }
    } else {
        cerr<<argv[0]<<" <INPUT_IMAGE>"<<endl;
        cerr<<argv[0]<<" <PRIMES_TO> <OUTPUT_FILE>"<<endl;
        return 1;
    }

    return 0;
}

