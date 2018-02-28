#pragma once

#include <Windows.h>

class Timer {
private:
	double dff;
	__int64 c1, c2;

public:
	Timer() {}

	void begin() {
		LARGE_INTEGER tmp;
		QueryPerformanceFrequency(&tmp);
		dff = tmp.QuadPart;
		QueryPerformanceCounter(&tmp);
		c1 = tmp.QuadPart;
	}

	// return ms
	double getDuration() {
		LARGE_INTEGER tmp;
		QueryPerformanceCounter(&tmp);
		c2 = tmp.QuadPart;
		return (c2 - c1) * 1000.0 / dff;
	}
};