/*
 * Debugging functions for GIES
 *
 * @author Alain Hauser
 * $Id: gies_debug.hpp 24 2012-01-05 13:51:28Z alhauser $
 */

#ifndef GIES_DEBUG_HPP_
#define GIES_DEBUG_HPP_

#include <iostream>

// Define default debug level
#ifndef DEBUG_OUTPUT_LEVEL
#define DEBUG_OUTPUT_LEVEL 0
#endif

namespace std {

/**
 * Sends a vector to an output stream (e.g., cout)
 */
template <typename T> ostream& operator<<(ostream& out, const vector<T>& vec) {
	int i;
	out << "(";
	for (i = 0; i < vec.size() - 1; i++)
		out << vec[i] << ", ";
	if (!vec.empty())
		out << vec.back();
	out << ")";
	return out;
}

/**
 * Sends a set to an output stream (e.g., cout)
 */
template <typename T> ostream& operator<<(ostream& out, const set<T>& s) {
	typename set<T>::iterator iter;
	out << "{";
	iter = s.begin();
	if (iter != s.end()) {
		out << *iter;
		for (++iter; iter != s.end(); ++iter)
			out << ", " << *iter;
	}
	out << "}";
	return out;
}

} // NAMESPACE std

// Macro for outputting debug messages depending on debug level
#if DEBUG_OUTPUT_LEVEL >= 1
#define DBOUT( level, message ) if ( DEBUG_OUTPUT_LEVEL >= level ) std::cout << message << std::endl
#else
#define DBOUT( level, message )
#endif

namespace std {

class nullbuf: public streambuf
{
	int xputc(int) { return 0; }
	streamsize xsputn(char const *, streamsize n) { return n; }
	int sync() { return 0; }
};

} // NAMESPACE std

class DebugStream
{
	int _level;
	std::nullbuf _nullbuf;
	std::ostream _nullstream;

public:
	DebugStream() : _level(0), _nullstream(&_nullbuf) {}
	DebugStream(int level) : _level(level), _nullstream(&_nullbuf) {}

	void setLevel(const int level) { _level = level; }

	int getLevel() const { return _level; }

	std::ostream& level(const int messageLevel)
	{
		if (messageLevel <= _level)
			return std::cout;
		return _nullstream;
	}
};

#ifdef DEFINE_GLOBAL_DEBUG_STREAM
#define GLOBAL
#else
#define GLOBAL extern
#endif

GLOBAL DebugStream dout;

#endif /* GIES_DEBUG_HPP_ */
