//************************************************************************
// --- CFEM++ Library
// ---
// --- Copyright 2024-2026 Ismaël Tchinda Ngueyong et al.
// --- ALL RIGHTS RESERVED
// ---
// --- This software is protected by international copyright laws.
// --- Use, distribution, or modification of this software in any form,
// --- source or binary, for personal, academic, or non-commercial 
// --- purposes is permitted free of charge, provided that this 
// --- copyright notice and this permission notice appear in all 
// --- copies and supporting documentation.
// ---
// --- The authors and contributors provide this code "as is" without 
// --- any express or implied warranty. In no event shall they be 
// --- held liable for any damages arising from the use of this software.
//************************************************************************

#ifndef PROFILER_HPP
#define PROFILER_HPP

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "ContractException.h"

namespace cfem
{
    struct ProfileData
    {
        std::string name;
        double durationMs;
        long memDiffKB;
        long totalMemKB;
    };

    class Profiler
    {
    private:
        std::vector<ProfileData> records;
        Profiler() {} // Private constructor for Singleton

    public:
        static Profiler &instance()
        {
            static Profiler inst;
            return inst;
        }

        // Moved to PUBLIC so ProfilingScope can see it
        size_t getRSS()
        {
            long rss = 0;
            std::ifstream fp("/proc/self/statm");
            if (fp.is_open())
            {
                fp >> rss; // Skip VmSize
                fp >> rss; // Read VmRSS (2nd field)
                return (size_t)rss * sysconf(_SC_PAGESIZE);
            }
            return 0;
        }

        // void addRecord(std::string name, double time, long mem)
        // {
        //     records.push_back({name, time, mem});
        // }

        void addRecord(std::string name, double time, long diff, long total) {
            records.push_back({name, time, diff, total});
        }

        void saveProfiler(std::string filename = "profile_results.csv") 
        {
            std::ofstream file(filename);
            if (!file.is_open()) return;

            file << "Method,TotalTime(ms),MemDiff(KB),TotalPeakMem(KB)\n";
            for (const auto &r : records) {
                file << r.name << "," << r.durationMs << "," << r.memDiffKB << "," << r.totalMemKB << "\n";
            }
            file.close();
#ifndef NDEBUG
            CFEM_INFO << "Profiler Exported to " << filename << ".\n";
#endif
        }
    };

    class ProfilingScope
    {
        std::string name;
        std::chrono::high_resolution_clock::time_point start;
        size_t startMem;

    public:
        ProfilingScope(std::string label) : name(label)
        {
            // Fix: Call getRSS directly from the instance
            startMem = Profiler::instance().getRSS();
            start = std::chrono::high_resolution_clock::now();
        }

        ~ProfilingScope() {
            auto end = std::chrono::high_resolution_clock::now();
            size_t endMem = Profiler::instance().getRSS();

            double duration = std::chrono::duration<double, std::milli>(end - start).count();
            long memDiff = (long)(endMem - startMem) / 1024;
            long totalPeak = (long)endMem / 1024; // Absolute memory usage

            Profiler::instance().addRecord(name, duration, memDiff, totalPeak);
        }
    };

// This helper forces the evaluation of __LINE__
#define PROFILE_SCOPE_CONCAT(name, line) ProfilingScope scope_##line(name)
#define PROFILE_SCOPE_HELPER(name, line) PROFILE_SCOPE_CONCAT(name, line)

// This is the macro you actually call
#define PROFILE_SCOPE(name) PROFILE_SCOPE_HELPER(name, __LINE__)

#endif
}
