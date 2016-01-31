#pragma once

#include <string>

using std::string;

enum LogLevel {
    LOG_INFO,
    LOG_OK,
    LOG_WARNING,
    LOG_ERROR
};

void logger(LogLevel level, string message);