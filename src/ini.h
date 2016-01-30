// Read an INI file into easy-to-access name/value pairs.

// inih and INIReader are released under the New BSD license (see LICENSE.txt).
// Go to the project home page for more info:
//
// http://code.google.com/p/inih/

/* Nonzero to allow multi-line value parsing, in the style of Python's
   ConfigParser. If allowed, ini_parse() will call the handler with the same
   name for each subsequent line parsed. */
#ifndef INI_ALLOW_MULTILINE
#define INI_ALLOW_MULTILINE 1
#endif

/* Nonzero to allow a UTF-8 BOM sequence (0xEF 0xBB 0xBF) at the start of
   the file. See http://code.google.com/p/inih/issues/detail?id=21 */
#ifndef INI_ALLOW_BOM
#define INI_ALLOW_BOM 1
#endif

/* Nonzero to use stack, zero to use heap (malloc/free). */
#ifndef INI_USE_STACK
#define INI_USE_STACK 1
#endif

/* Stop parsing on first error (default is to keep parsing). */
#ifndef INI_STOP_ON_FIRST_ERROR
#define INI_STOP_ON_FIRST_ERROR 0
#endif

/* Maximum line length for any line in INI file. */
#ifndef INI_MAX_LINE
#define INI_MAX_LINE 200
#endif

#include <map>
#include <set>
#include <string>

// Read an INI file into easy-to-access name/value pairs. (Note that I've gone
// for simplicity here rather than speed, but it should be pretty decent.)
class INIReader {
  public:
    // Construct INIReader and parse given filename. See ini.h for more info
    // about the parsing.
    INIReader(std::string filename);
    ~INIReader();

    // Return the result of ini_parse(), i.e., 0 on success, line number of
    // first error on parse error, or -1 on file open error.
    int ParseError();

    // Get a string value from INI file, returning default_value if not found.
    std::string Get(std::string section, std::string name,
                    std::string default_value);

    // Get an integer (long) value from INI file, returning default_value if
    // not found or not a valid integer (decimal "1234", "-1234", or hex
    // "0x4d2").
    long GetInteger(std::string section, std::string name, long default_value);

    // Get a real (floating point double) value from INI file, returning
    // default_value if not found or not a valid floating point value
    // according to strtod().
    double GetReal(std::string section, std::string name, double default_value);

    // Get a boolean value from INI file, returning default_value if not found
    // or if
    // not a valid true/false value. Valid true values are "true", "yes", "on",
    // "1",
    // and valid false values are "false", "no", "off", "0" (not case
    // sensitive).
    bool GetBoolean(std::string section, std::string name, bool default_value);

    // Returns all the section names from the INI file, in alphabetical order,
    // but in the
    // original casing
    std::set<std::string> GetSections() const;

    // Returns all the field names from a section in the INI file, in
    // alphabetical order,
    // but in the original casing. Returns an empty set if the field name is
    // unknown
    std::set<std::string> GetFields(std::string section) const;

  private:
    int _error;
    std::map<std::string, std::string> _values;
    // Because we want to retain the original casing in _fields, but
    // want lookups to be case-insensitive, we need both _fields and _values
    std::set<std::string> _sections;
    std::map<std::string, std::set<std::string> *> _fields;
    static std::string MakeKey(std::string section, std::string name);
    static int ValueHandler(void * user, const char * section,
                            const char * name, const char * value);
};
