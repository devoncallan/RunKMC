#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <sstream>

#include "core/C.h"
#include "utils/parse.h"

namespace console
{
    namespace color = C::io::color;

    enum class LogLevel : uint8_t
    {
        SILENT = 0,
        ERROR = 1,
        WARNING = 2,
        INFO = 3,
        DEBUG = 4,
        TRACE = 5
    };

    class NullStream : public std::ostream
    {
        class NullBuffer : public std::streambuf
        {
        public:
            int overflow(int c) override { return c; }
        } buffer;

    public:
        NullStream() : std::ostream(&buffer) {}
    };

    class Logger
    {
    private:
        static inline LogLevel currentLevel = LogLevel::INFO;
        static inline NullStream nullStream;
        static inline std::ofstream logFile;
        static inline bool fileLoggingEnabled = false;

    public:
        static void setLevel(LogLevel level) { currentLevel = level; }

        static void enableFileLogging(const std::filesystem::path &path)
        {
            logFile.open(path, std::ios::out | std::ios::app);
            fileLoggingEnabled = logFile.is_open();
            if (!fileLoggingEnabled)
                std::cout << "Failed to open log file: " << path.string() << std::endl;
            // warning("Failed to open log file: " + path.string());
        }

        static void disableFileLogging()
        {
            if (logFile.is_open())
                logFile.close();
            fileLoggingEnabled = false;
        }

        static bool shouldLog(LogLevel level) { return level <= currentLevel; }
        static std::ostream &consoleStream() { return std::cout; }
        static std::ostream &fileStream() { return logFile; }
        static bool isFileLoggingEnabled() { return fileLoggingEnabled; }
    };

    static inline void _printWithContext(
        std::string_view title,
        std::string_view context,
        std::string_view msg,
        std::string_view c,
        LogLevel level,
        bool fatal = false)
    {
        if (!Logger::shouldLog(level))
            return;

        std::cout << color::on(c);
        std::cout << "[" << context << "] ";
        std::cout << "[" << title << "] : ";
        std::cout << msg;
        std::cout << color::on(color::DEFAULT) << std::endl;

        if (Logger::isFileLoggingEnabled())
        {
            auto &fs = Logger::fileStream();
            fs << "[" << context << "] ";
            fs << "[" << title << "] : ";
            fs << msg << std::endl;
            fs.flush();
        }

        if (fatal)
            exit(EXIT_FAILURE);
    }

    class LogContext
    {
    private:
        const char *name;

    public:
        constexpr LogContext(const char *n) : name(n) {}
        inline void trace(std::string_view msg) const
        {
            if (!Logger::shouldLog(LogLevel::TRACE))
                return;
            _printWithContext("TRACE", name, msg, color::CYN, LogLevel::TRACE);
        }

        inline void debug(std::string_view msg) const
        {
            if (!Logger::shouldLog(LogLevel::DEBUG))
                return;
            _printWithContext("DEBUG", name, msg, color::BLU, LogLevel::DEBUG);
        }

        inline void info(std::string_view msg) const
        {
            if (!Logger::shouldLog(LogLevel::INFO))
                return;
            _printWithContext("INFO", name, msg, color::GRN, LogLevel::INFO);
        }

        inline void warning(std::string_view msg) const
        {
            if (!Logger::shouldLog(LogLevel::WARNING))
                return;
            _printWithContext("WARNING", name, msg, color::YLW, LogLevel::WARNING);
        }

        [[noreturn]] void error(std::string_view msg) const
        {
            _printWithContext("ERROR", name, msg, color::RED, LogLevel::ERROR, true);
        }
    };
}

namespace console::term
{

    inline constexpr std::string_view ESC = "\x1b";
    inline constexpr std::string_view ST = "\x1b\\";                 // String Terminator
    inline constexpr std::string_view OSC = "\x1b]";                 // Operating System Command
    inline constexpr std::string_view OSC8_OPEN = "\x1b]8;;";        // OSC 8 "open link"
    inline constexpr std::string_view OSC8_CLOSE = "\x1b]8;;\x1b\\"; // OSC 8 "close link"

    // Write an OSC 8 hyperlink to the given output stream
    inline std::ostream &_osc8(std::ostream &out, std::string_view uri, std::string_view label)
    {
        out << OSC8_OPEN << uri << ST;

        if (!label.empty())
            out << label;
        out << OSC8_CLOSE;
        return out;
    }

    // Return an OSC 8 hyperlink string
    inline std::string osc8(std::string_view uri, std::string_view label)
    {
        std::ostringstream oss;
        _osc8(oss, uri, label);
        return oss.str();
    }

    // Return an OSC 8 hyperlink string for a filesystem path
    inline std::string linkPath(const std::filesystem::path &p, std::string_view label = {})
    {
        // POSIX file:// URI
        const std::string fileUri = std::string("file://") + p.string();
        const std::string labelStr = label.empty() ? p.filename().string() : std::string(label);
        return osc8(fileUri, labelStr);
    }
};

namespace console
{
    // General console message with optional color and fatal flag
    static inline void _print(std::string_view title, std::string_view msg, std::string_view c = color::DEFAULT, bool fatal = false)
    {
        std::cout << color::on(c) << "[" << title << "] : " << msg << color::on(color::DEFAULT) << std::endl;
        if (fatal)
            exit(EXIT_FAILURE);
    }

    static inline void _printVector(std::string_view title, const std::vector<std::string_view> &vec, std::string_view c = color::DEFAULT)
    {
        std::cout << color::on(c) << "[" << title << "] : " << str::join(vec, ", ") << color::on(color::DEFAULT) << std::endl;
    }

    static inline void _printLink(std::string_view title, const std::filesystem::path &path, std::string_view label = {}, std::string_view c = color::DEFAULT)
    {
        _print(title, term::linkPath(path, label), c);
    }

    static inline void debug(std::string_view msg) { _print("DEBUG", msg, color::BLU); }
    static inline void debug_msg(std::string_view msg) { _print("DEBUG", msg, color::BLU); }
    static inline void log(std::string_view msg) { _print("LOG", msg, color::GRN); }
    static inline void warning(std::string_view msg) { _print("WARNING", msg, color::YLW); }
    static inline void input_warning(std::string_view msg) { _print("INPUT WARNING", msg, color::YLW); }
    static inline void error(std::string_view msg) { _print("ERROR", msg, color::RED, true); }
    static inline void input_error(std::string_view msg) { _print("INPUT ERROR", msg, color::RED, true); }
    static inline void link(const std::filesystem::path &path, std::string_view label = {}) { _printLink("LINK", path, label, color::GRN); }

    static inline void logVector(const std::vector<std::string_view> &vec) { log(str::join(vec, ", ")); }
    static inline void logVector(const std::vector<std::string> &vec) { log(str::join(vec, ", ")); }
};
