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
#include <stdexcept>

namespace cfem {

/**
 * @brief Custom exception for CFEM errors
 */
class CfemException : public std::runtime_error {
public:
    explicit CfemException(const std::string& msg) : std::runtime_error(msg) {}
};

enum class LogLevel { Info, Warning, Error };

class Logger {
public:
    static Logger& getInstance() {
        static Logger instance;
        return instance;
    }

    struct LogStream {
        std::stringstream ss;
        LogLevel level;
        const char* file;
        int line;

        LogStream(LogLevel l, const char* f, int ln) : level(l), file(f), line(ln) {}

        /**
         * @brief Destructor handles printing and throwing exceptions for Error level
         */
        ~LogStream() noexcept(false) { // Note: noexcept(false) allows throwing from destructor
            std::string prefix;
            std::string color;
            bool showMetadata = true;

            switch (level) {
                case LogLevel::Info:
                    prefix = "[INFO] ";
                    color = "\033[1;32m"; // Green
                    showMetadata = false;
                    break;
                case LogLevel::Warning:
                    prefix = "[WARNING] ";
                    color = "\033[1;33m"; // Yellow
                    break;
                case LogLevel::Error:
                    prefix = "[ERROR] ";
                    color = "\033[1;31m"; // Red
                    break;
            }
            
            // Output the message
            std::cerr << color << prefix << "\033[0m" << ss.str();
            
            if (showMetadata) {
                std::cerr << " (" << file << ":" << line << ")";
            }
            
            std::cerr << std::endl;

            // --- The Critical Change ---
            if (level == LogLevel::Error) {
                throw CfemException(ss.str());
            }
        }

        template <typename T>
        LogStream& operator<<(const T& msg) {
            ss << msg;
            return *this;
        }
    };

    LogStream build(LogLevel level, const char* file, int line) {
        return LogStream(level, file, line);
    }

private:
    Logger() = default;
};

} // namespace cfem

#define CFEM_INFO  cfem::Logger::getInstance().build(cfem::LogLevel::Info, __FILE__, __LINE__)
#define CFEM_WARN  cfem::Logger::getInstance().build(cfem::LogLevel::Warning, __FILE__, __LINE__)
#define CFEM_ERROR cfem::Logger::getInstance().build(cfem::LogLevel::Error, __FILE__, __LINE__)