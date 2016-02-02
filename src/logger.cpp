#include <iostream>
#include "logger.h"

using std::cout;
using std::flush;

void logger(LogLevel level, string message) {
    string set_color, reset_color = "\033[0;0m";
    switch (level) {
    case LOG_OK:
        set_color = "\033[0;32m";
        break;
    case LOG_WARNING:
        set_color = "\033[0;33m";
        break;
    case LOG_ERROR:
        set_color = "\033[0;31m";
        break;
    default:
        set_color = reset_color;
    }
    cout << set_color << message << reset_color << flush;
}