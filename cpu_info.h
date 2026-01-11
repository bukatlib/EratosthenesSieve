#ifndef HLIDAC_PER_CPU_INFO_H
#define HLIDAC_PER_CPU_INFO_H

#include <filesystem>
#include <set>
#include <string>
#include <optional>
#include <stdint.h>
#include <vector>

enum class CpuArch { AArch64, AMD64 };
enum class CacheType { Instruction, Data, Unified };
enum class Aggregate { Min, Max };

std::ostream& operator<<(std::ostream& os, const CpuArch& ca);
std::ostream& operator<<(std::ostream& os, const CacheType& ca);

using prop_ca = std::optional<CpuArch>;
using prop_ct = std::optional<CacheType>;
using prop_ui = std::optional<uint64_t>;

struct CacheInfo {
    public:
        CacheInfo(std::filesystem::path cache_dir);
        bool operator<(const CacheInfo& ci) const;
        void show() const;

        prop_ui level() const;
        prop_ui size() const;
        prop_ui number_of_sets() const;
        prop_ui ways_of_associativity() const;
        prop_ui coherency_line_size() const;
        prop_ui shared_cpu_count() const;
        prop_ct type() const;

    private:
        template <class T>
        static void show_helper(const std::string& label, const std::optional<T>& opt);

        prop_ui level_p;
        prop_ui size_p;
        prop_ui number_of_sets_p;
        prop_ui ways_of_associativity_p;
        prop_ui coherency_line_size_p;
        prop_ui shared_cpu_count_p;
        prop_ct type_p;

        static std::set<std::string> check_files;
};

struct CoreInfo {
    public:
        CoreInfo(uint64_t core_idx, std::filesystem::path cpu_dir);
        void show() const;
        std::vector<CacheInfo> caches() const;

    private:
        uint64_t core_idx;
        std::vector<CacheInfo> caches_p;
};


struct CpuInfo {
    public:
        CpuInfo();
        void show() const;
        std::vector<CoreInfo> cores() const;
        prop_ca arch() const;

        prop_ui cache_size(uint64_t level, CacheType type, Aggregate fce = Aggregate::Min) const;
        prop_ui cache_size_exclusive(uint64_t level, CacheType type, Aggregate fce = Aggregate::Min) const;

    private:
        static prop_ca setup_arch();
        prop_ui cache_size_helper(uint64_t level, CacheType type, Aggregate fce, bool exclusive) const;

        prop_ca arch_p;    
        std::vector<CoreInfo> cores_p;
};

#endif
