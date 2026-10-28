/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2026                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef VCG_WRAP_SYSTEM_UTF8_FILE_H
#define VCG_WRAP_SYSTEM_UTF8_FILE_H

/** UTF-8 aware file access.

  Every filename taken by the VCG importers/exporters is a UTF-8 encoded
  const char*. On POSIX systems the bytes are passed through unchanged; on
  Windows they are converted to UTF-16 and the wide CRT/stream entry points are
  used, since the narrow ones interpret paths in the active ANSI code page.

  For backward compatibility, a Windows path that is not valid UTF-8 is assumed
  to be in the ANSI code page and is opened with the narrow API, as before.
*/

#include <cstdio>
#include <ctime>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#include <filesystem>
#include <io.h>
#else
#include <unistd.h>
#endif

namespace vcg {
namespace utf8 {

/// Decodes a NUL-terminated UTF-8 string into UTF-16 code units.
/// Returns false (leaving out unspecified) on malformed input: truncated or
/// overlong sequences, surrogate code points, values above U+10FFFF.
inline bool ToUtf16(const char *s, std::wstring &out)
{
	out.clear();
	const unsigned char *p = reinterpret_cast<const unsigned char *>(s);
	while (*p) {
		unsigned int cp;
		int extra;
		unsigned int c = *p++;
		if (c < 0x80)      { cp = c;        extra = 0; }
		else if (c < 0xC2) return false; // continuation byte or overlong 2-byte lead
		else if (c < 0xE0) { cp = c & 0x1F; extra = 1; }
		else if (c < 0xF0) { cp = c & 0x0F; extra = 2; }
		else if (c < 0xF5) { cp = c & 0x07; extra = 3; }
		else return false;
		for (int i = 0; i < extra; ++i) {
			unsigned int cc = *p;
			if ((cc & 0xC0) != 0x80) return false; // also catches the terminator
			cp = (cp << 6) | (cc & 0x3F);
			++p;
		}
		if ((extra == 2 && cp < 0x800) || (extra == 3 && cp < 0x10000) ||
		    cp > 0x10FFFF || (cp >= 0xD800 && cp <= 0xDFFF))
			return false;
		if (cp < 0x10000)
			out.push_back(wchar_t(cp));
		else {
			cp -= 0x10000;
			out.push_back(wchar_t(0xD800 + (cp >> 10)));
			out.push_back(wchar_t(0xDC00 + (cp & 0x3FF)));
		}
	}
	return true;
}

#ifdef _WIN32
/// Path type accepted by std::fstream constructors and open().
typedef std::filesystem::path StreamPath;
#else
typedef std::string StreamPath;
#endif

/// Converts a UTF-8 filename to something std::ifstream / std::ofstream can open:
///   std::ifstream in(vcg::utf8::ToStreamPath(filename));
inline StreamPath ToStreamPath(const char *filename)
{
#ifdef _WIN32
	std::wstring w;
	if (ToUtf16(filename, w))
		return std::filesystem::path(w);
	return std::filesystem::path(filename); // legacy ANSI path
#else
	return std::string(filename);
#endif
}

/// fopen() replacement taking a UTF-8 filename.
inline FILE *FOpen(const char *filename, const char *mode)
{
#ifdef _WIN32
	std::wstring wname, wmode;
	if (ToUtf16(filename, wname) && ToUtf16(mode, wmode))
		return _wfopen(wname.c_str(), wmode.c_str());
#endif
	return std::fopen(filename, mode);
}

/// access() replacement taking a UTF-8 filename. Mode follows access()/_access()
/// (0 existence, 2 write, 4 read, 6 read/write).
inline int Access(const char *filename, int mode)
{
#ifdef _WIN32
	std::wstring w;
	if (ToUtf16(filename, w))
		return _waccess(w.c_str(), mode);
	return _access(filename, mode);
#else
	return access(filename, mode);
#endif
}

/// mkdir() replacement taking a UTF-8 path. Returns 0 on success.
inline int MkDir(const char *dirname)
{
#ifdef _WIN32
	std::wstring w;
	if (ToUtf16(dirname, w))
		return _wmkdir(w.c_str());
	return _mkdir(dirname);
#else
	return mkdir(dirname, 0755);
#endif
}

/// Reads the last modification time of a file given its UTF-8 name.
inline bool ModificationTime(const char *filename, time_t &mtime)
{
#ifdef _WIN32
	std::wstring w;
	if (ToUtf16(filename, w)) {
		struct _stat64 st;
		if (_wstat64(w.c_str(), &st) != 0) return false;
		mtime = st.st_mtime;
		return true;
	}
	struct _stat64 st;
	if (_stat64(filename, &st) != 0) return false;
	mtime = st.st_mtime;
	return true;
#else
	struct stat st;
	if (stat(filename, &st) != 0) return false;
	mtime = st.st_mtime;
	return true;
#endif
}

} // namespace utf8
} // namespace vcg

#endif // VCG_WRAP_SYSTEM_UTF8_FILE_H
