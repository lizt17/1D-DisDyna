#pragma once
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <cctype>

// 去掉字符串首尾空白字符
inline std::string trim(const std::string& s)
{
    size_t begin = 0;
    while (begin < s.size() && std::isspace(static_cast<unsigned char>(s[begin])))
    {
        ++begin;
    }

    size_t end = s.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(s[end - 1])))
    {
        --end;
    }

    return s.substr(begin, end - begin);
}

// 去掉注释，支持：
// 1) # comment
// 2) // comment
inline std::string removeComment(const std::string& line)
{
    size_t posHash  = line.find('#');
    size_t posSlash = line.find("//");

    size_t cutPos = std::string::npos;

    if (posHash != std::string::npos)
    {
        cutPos = posHash;
    }

    if (posSlash != std::string::npos)
    {
        if (cutPos == std::string::npos)
        {
            cutPos = posSlash;
        }
        else
        {
            cutPos = std::min(cutPos, posSlash);
        }
    }

    if (cutPos != std::string::npos)
    {
        return line.substr(0, cutPos);
    }

    return line;
}

// 从 input.txt 中读取形如 key=value; 的参数
// 支持注释和空行
inline std::map<std::string, double> readFile(const std::string& filename)
{
    std::map<std::string, double> params;

    std::ifstream infile(filename.c_str());
    if (!infile)
    {
        std::cerr << "Error: Cannot open input file: " << filename << std::endl;
        return params;
    }

    std::string line;
    int lineNo = 0;

    while (std::getline(infile, line))
    {
        ++lineNo;

        // 去掉注释
        line = removeComment(line);

        // 去掉首尾空格
        line = trim(line);

        // 跳过空行
        if (line.empty())
        {
            continue;
        }

        std::stringstream ss(line);
        std::string token;

        // 支持一行多个 key=value; 对
        while (std::getline(ss, token, ';'))
        {
            token = trim(token);
            if (token.empty())
            {
                continue;
            }

            size_t eq_pos = token.find('=');
            if (eq_pos == std::string::npos)
            {
                std::cerr << "Warning: line " << lineNo
                          << " ignored, missing '=' : " << token << std::endl;
                continue;
            }

            std::string key    = trim(token.substr(0, eq_pos));
            std::string valStr = trim(token.substr(eq_pos + 1));

            if (key.empty() || valStr.empty())
            {
                std::cerr << "Warning: line " << lineNo
                          << " ignored, empty key/value : " << token << std::endl;
                continue;
            }

            char* endPtr = NULL;
            double value = std::strtod(valStr.c_str(), &endPtr);

            if (endPtr == valStr.c_str() || *endPtr != '\0')
            {
                std::cerr << "Warning: line " << lineNo
                          << " ignored, invalid number for key [" << key
                          << "] : " << valStr << std::endl;
                continue;
            }

            params[key] = value;
        }
    }

    return params;
}

// 将参数写入 outputVars.txt 文件，格式为 key = value
inline void writeFile(const std::string& filename, const std::map<std::string, double>& params)
{
    std::ofstream fout(filename.c_str());
    if (!fout)
    {
        std::cerr << "Error: Cannot write to file: " << filename << std::endl;
        return;
    }

    for (std::map<std::string, double>::const_iterator it = params.begin(); it != params.end(); ++it)
    {
        fout << it->first << " = " << it->second << std::endl;
    }
}