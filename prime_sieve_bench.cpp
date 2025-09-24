#include <chrono>
#include <iostream>
#include <string>
#include <stdint.h>
#include "eratosthenes_sieve.h"

using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[])  {
	uint64_t to_value;
    if (argc >= 2)  {
        string to_value_str = argv[1];
        string::size_type sz = 0; 
        to_value = stoull(to_value_str, &sz);
        if (sz != to_value_str.size())  {
            cerr<<"Invalid '<PRIMES_TO>' argument!"<<endl;
            cerr<<argv[0]<<" <PRIMES_TO> <OUTPUT_FILE>"<<endl;
            return 1;
        }
    } else {
        bool valid;
        do {
            valid = true;
            string to_value_str;

            cin.clear();
	        cout<<"Primes to: ";
            cin>>to_value_str;

            if (cin.eof())  {
                cerr<<"No answer provided, terminating..."<<endl;  
                return 1;
            }

            if (to_value_str.empty() || to_value_str[0] == '-')
                valid = false;

            try {
                string::size_type sz = 0; 
                to_value = stoull(to_value_str, &sz);
                if (sz != to_value_str.size())
                    valid = false;
            } catch (...)   {
                valid = false;
            }

            if (!valid)
                cerr<<"Invalid input, try again!"<<endl;

	    } while (!valid);

        // Ignore '\n' character at the end.
        cin.clear();
        cin.ignore(1, '\n');
    }

	string file_name;
    if (argc >= 3)  {
        string arg_output = argv[2];
        if (arg_output != "-" && !arg_output.empty())
            file_name = arg_output;
    } else {
	    cout<<"Output file (empty no write): ";
        getline(cin, file_name);
	    cout<<string(60,'#')<<endl;
    }


    time_point<system_clock, duration<double>> tp1 = high_resolution_clock::now();

    Eratosthenes sieve(to_value);
    cout<<"sieve size: "<<sieve.bit_size()<<" bits ("<<sieve.bytes_allocated()<<" bytes)"<<endl;

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
        sieve.write_primes_to_file(file_name);
        time_point<system_clock, duration<double>> tp5 = high_resolution_clock::now();
        double write_time = duration_cast<duration<double>>(tp5 - tp4).count();
        cout<<"Primes written to '"<<file_name<<"' file in "<<write_time<<" secs."<<endl;
    }

    return 0;
}

