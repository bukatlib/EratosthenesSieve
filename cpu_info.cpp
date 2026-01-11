#include <algorithm>
#include <fstream>
#include <iostream>
#include "cpu_info.h"
#include "cmake_vars.h"

using namespace std;
using namespace std::filesystem;

ostream& operator<<(ostream& os, const CpuArch& ca) {
    switch (ca) {
        case CpuArch::AArch64:
            os<<"aarch64";
            break;
        case CpuArch::AMD64:
            os<<"amd64";
            break;
    }
    return os;
}

ostream& operator<<(ostream& os, const CacheType& ct)   {
    switch (ct) {
        case CacheType::Instruction:
            os<<"instruction";
            break;
        case CacheType::Data:
            os<<"data";
            break;
        case CacheType::Unified:
            os<<"unified";
            break;
    }
    return os;
}

CacheInfo::CacheInfo(path cache_dir)  {
    for (const directory_entry& entry : directory_iterator(cache_dir))   {
        if (!entry.is_regular_file())
            continue;

        path file_path = entry.path();
        string filename = file_path.filename().string();

        if (!CacheInfo::check_files.contains(filename))
            continue;

        string record;
        ifstream cache_file(file_path.string(), ios::in);
        if (!cache_file)    {
            cerr<<"Failed to read "<<entry<<" file, skipping!"<<endl;
            continue;
        }
        cache_file>>record;
        cache_file.close();
        
        if (record.empty()) {
            cerr<<"Read empty record for "<<filename<<", skipping..."<<endl;
            continue;
        }

        if (filename == "level") {
            level_p = stoull(record);
        } else if (filename == "size")  {
            char unit = record.back();
            record.pop_back();
            uint64_t cache_size = stoull(record);
            switch (unit)   {
                case 'K':
                    size_p = 1024ul * cache_size;
                    break;
                case 'M':
                    size_p = 1024ul * 1024ul * cache_size;
                    break;
                default:
                    cerr<<"Unit not detected in "<<filename<<", skipping..."<<endl;
                    continue;
            }
        } else if (filename == "number_of_sets")    {
            number_of_sets_p = stoull(record);
        } else if (filename == "ways_of_associativity")  {
            ways_of_associativity_p = stoull(record);
        } else if (filename == "coherency_line_size")   {
            coherency_line_size_p = stoull(record);
        } else if (filename == "shared_cpu_map")    {
            uint64_t shared_cpu_map = stoul("0x" + record, nullptr, 16);
            shared_cpu_count_p = popcount(shared_cpu_map);
        } else if (filename == "type")  {
            if (record == "Instruction")
                type_p = CacheType::Instruction;
            else if (record == "Data")
                type_p = CacheType::Data;
            else if (record == "Unified")
                type_p = CacheType::Unified;
            else
                cerr<<"Unknown cache type '"<<record<<"'!"<<endl;
        }
    }
}

bool CacheInfo::operator<(const CacheInfo& ci) const    {
    if (level_p < ci.level_p)
        return true;
    else if (level_p > ci.level_p)
        return false;
    else
        return type_p < ci.type_p;
}

void CacheInfo::show() const    {
    CacheInfo::show_helper("Level", level_p);
    CacheInfo::show_helper("Type", type_p);
    CacheInfo::show_helper("Size", size_p);
    CacheInfo::show_helper("Number of sets", number_of_sets_p);
    CacheInfo::show_helper("Ways of associativity", ways_of_associativity_p);
    CacheInfo::show_helper("Coherency line size", coherency_line_size_p);
    CacheInfo::show_helper("Shared count", shared_cpu_count_p);
}

prop_ui CacheInfo::level() const    {
    return level_p;
}

prop_ui CacheInfo::size() const {
    return size_p;
}

prop_ui CacheInfo::number_of_sets() const   {
    return number_of_sets_p;
}

prop_ui CacheInfo::ways_of_associativity() const    {
    return ways_of_associativity_p;
}

prop_ui CacheInfo::coherency_line_size() const  {
    return coherency_line_size_p;
}

prop_ui CacheInfo::shared_cpu_count() const {
    return shared_cpu_count_p;
}

prop_ct CacheInfo::type() const {
    return type_p;
}

template <class T>
void CacheInfo::show_helper(const std::string& label, const std::optional<T>& opt)    {
    if (opt)
        cout<<label<<": "<<opt.value()<<endl;
}

set<string> CacheInfo::check_files = {
    "level", "size", "number_of_sets", "ways_of_associativity", "coherency_line_size", "shared_cpu_map", "type"
};

CoreInfo::CoreInfo(uint64_t core_idx, path cpu_dir) : core_idx(core_idx) {
    uint64_t cache_idx = 0ul;
    path cache_prefix = cpu_dir / "cache";
    path cache_dir = cache_prefix / string("index" + to_string(cache_idx));
    while (exists(cache_dir) && is_directory(cache_dir))    {
        caches_p.emplace_back(cache_dir);
        cache_dir = cache_prefix / string("index" + to_string(++cache_idx));
    }

    sort(caches_p.begin(), caches_p.end());
}

void CoreInfo::show() const    {
    cout<<"core idx: "<<core_idx<<endl<<endl;
    for (const CacheInfo& ci : caches_p)  {
        ci.show();
        cout<<endl;
    }
}

vector<CacheInfo> CoreInfo::caches() const  {
    return caches_p;
}

CpuInfo::CpuInfo() : arch_p(setup_arch()) {
    uint64_t core_idx = 0ul;
    string cpu_path_prefix = "/sys/devices/system/cpu/cpu";
    path cpu_dir_path = cpu_path_prefix + to_string(core_idx);
    while (exists(cpu_dir_path) && is_directory(cpu_dir_path))    {
        cores_p.emplace_back(core_idx, cpu_dir_path);
        cpu_dir_path = string(cpu_path_prefix) + to_string(++core_idx);
    }
}

void CpuInfo::show() const    {
    cout<<"Processor model: "<<PROCESSOR_NAME<<endl;
    if (arch_p)
        cout<<"CPU architecture: "<<arch_p.value()<<endl;
    cout<<"Threads per core: "<<THREADS_PER_CORE<<endl;
    cout<<"64 bit support: "<<IS_64BIT<<endl;
    cout<<"Physical memory: "<<TOTAL_PHYSICAL_MEMORY<<" MB"<<endl;
    cout<<"Operating system: "<<OPERATING_SYSTEM<<endl<<endl;
    
    for (const CoreInfo& ci : cores_p)    {
        cout<<string(30, '#')<<endl;
        ci.show();
    }
}

vector<CoreInfo> CpuInfo::cores() const {
    return cores_p;
}

prop_ca CpuInfo::arch() const {
    return arch_p;
}

prop_ui CpuInfo::cache_size(uint64_t level, CacheType type, Aggregate fce) const  {
    return cache_size_helper(level, type, fce, false);
}

prop_ui CpuInfo::cache_size_exclusive(uint64_t level, CacheType type, Aggregate fce) const  {
    return cache_size_helper(level, type, fce, true);
}

prop_ca CpuInfo::setup_arch() {
    #if defined(__aarch64__)
    return CpuArch::AArch64;
    #elif defined(__x86_64__)
    return CpuArch::AMD64;
    #else
    return nullopt;
    #endif
}

prop_ui CpuInfo::cache_size_helper(uint64_t level, CacheType type, Aggregate fce, bool exclusive) const  {
    vector<uint64_t> sizes;
    for (const CoreInfo& ci : cores_p)  {
        for (const CacheInfo& cache : ci.caches())  {
            bool level_has_size = (cache.level() == level) && cache.size();
            bool type_match = (cache.type() == type) || (cache.type() == CacheType::Unified);
            if (level_has_size && type_match)   {
                uint64_t cache_size = cache.size().value();
                if (exclusive && cache.shared_cpu_count())  {
                    uint64_t shared_count = cache.shared_cpu_count().value();
                    if (shared_count > 0ul)
                        cache_size /= shared_count;
                }
                sizes.push_back(cache_size);
            }
        }
    }

    if (!sizes.empty()) {
        if (fce == Aggregate::Min)
            return *min_element(sizes.cbegin(), sizes.cend());
        else
            return *max_element(sizes.cbegin(), sizes.cend());
    } else {
        return nullopt;
    }
}

