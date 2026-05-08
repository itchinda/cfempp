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

#pragma once
#include <iostream>
#include <sstream>
#include <string>
#include <atomic>
#include <vector>
#include "ContractException.h" // Intégration de tes contrats

namespace cfem {

/**
 * @brief État global des erreurs (Lock-Free / HPC-Ready)
 */
class ErrorManager {
public:
    static ErrorManager& instance() {
        static ErrorManager inst;
        return inst;
    }

    void addWarning() noexcept { warnings_.fetch_add(1, std::memory_order_relaxed); }
    void addError() noexcept   { errors_.fetch_add(1, std::memory_order_relaxed); }

    bool hasErrors() const noexcept { return errors_.load(std::memory_order_acquire) > 0; }
    uint32_t errorCount() const noexcept { return errors_.load(std::memory_order_relaxed); }

    void clear() noexcept {
        errors_.store(0);
        warnings_.store(0);
    }

private:
    ErrorManager() = default;
    std::atomic<uint32_t> errors_{0};
    std::atomic<uint32_t> warnings_{0};
};

enum class LogLevel { Info, Warning, Error, Fatal };

class Logger {
public:
    struct LogStream {
        std::stringstream ss;
        LogLevel level;
        const char* file;
        int line;

        LogStream(LogLevel l, const char* f, int ln) : level(l), file(f), line(ln) {}

        ~LogStream() {
            std::string msg = ss.str();
            const char* color = "";
            const char* prefix = "";

            switch (level) {
                case LogLevel::Info:    
                    prefix = "[INFO] ";    
                    color = "\033[32m"; 
                    break;
                case LogLevel::Warning: 
                    prefix = "[WARN] "; 
                    color = "\033[33m"; 
                    ErrorManager::instance().addWarning(); break;
                case LogLevel::Error:   
                    prefix = "[ERROR] ";
                    color = "\033[31m"; 
                    ErrorManager::instance().addError(); break;
                case LogLevel::Fatal:   
                    prefix = "[FATAL] ";
                    color = "\033[1;31m"; break;
            }

            std::cerr << color << prefix << "\033[0m" << msg;
#ifndef NDEBUG
            if (level != LogLevel::Info) 
                std::cerr << " (" << file << ":" << line << ")";
#endif
            std::cerr << std::endl;

            if (level == LogLevel::Fatal) std::abort(); 
        }

        template<typename T>
        LogStream& operator<<(const T& v) { ss << v; return *this; }
    };

    static LogStream build(LogLevel l, const char* f, int ln) { return LogStream(l, f, ln); }
};

} // namespace cfem

// --- Macros de Logging (Non-bloquant pour Error/Warning)
#define CFEM_INFO  cfem::Logger::build(cfem::LogLevel::Info, __FILE__, __LINE__)
#define CFEM_WARN  cfem::Logger::build(cfem::LogLevel::Warning, __FILE__, __LINE__)
#define CFEM_ERROR cfem::Logger::build(cfem::LogLevel::Error, __FILE__, __LINE__)
#define CFEM_FATAL cfem::Logger::build(cfem::LogLevel::Fatal, __FILE__, __LINE__)

#define CFEM_THROW(msg) \
    do { \
        std::ostringstream ss_; \
        ss_ << msg; \
        cfem::ErrorManager::instance().addError(); \
        std::string tmp_msg = ss_.str(); \
        throw cfem::ContractException(__FILE__, __LINE__, "CFEM_THROW", "Runtime Error", tmp_msg.c_str()); \
    } while(0)
