#ifndef REGIMPLICIT_CONSTANTS_H
#define REGIMPLICIT_CONSTANTS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <cassert>
#include <cmath>
#include <chrono>

using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

/* ===== Physical Constants ===== */
const float c0   = 299792458.0f;     // m/s
const float eps0 = 8.85418782e-12f;  // F/m
const float mu0  = 1.25663706e-6f;   // H/m
const float pi = 3.1415926535f;

const float freq   = 14.0e6f;
const float lambda = c0 / freq;
const float T      = 1.0f / freq;
const float t0     = 3.0f * T;
const float sig0   = 1.0f * T;
const float n0     = 2.0f;

const float ez0    = 1.0f; // V/m

/* ====== Spatial Constants ====== */
const int x_resolution = 16;
const int num_wavelengths = 6;

const int nx = num_wavelengths * x_resolution + 1;
const int ny = nx;
const int nz = nx;

const float dx = (num_wavelengths * lambda - 0.0f) / float(nx - 1);
const float dy = dx;
const float dz = dx;

/* ====== Temporal Constants ====== */
// Courant number
const float cfl = 0.35f;
const float dt = cfl * dx / c0;
const int nt = int(num_wavelengths * lambda / (c0 * dt)) + 1;

/* ====== Derivative Coefficients ====== */
const float ddx = 1.0f / dx;
const float ddy = 1.0f / dy;
const float ddz = 1.0f / dz;

const float c1 = dt / (2.0f * eps0);
const float c2 = dt / (2.0f * mu0);

/* ====== Helpful Constants ====== */
const int half_nx = int(nx / 2);
const int half_ny = int(ny / 2);
const int half_nz = int(ny / 2);

/* ===== Timer struct ===== */
struct Timer {
	steady_clock::time_point start_time;
	steady_clock::time_point end_time;
	steady_clock::time_point last;

	float elapsed = 0.0f;
	float total = 0.0f;

	void start() {
		start_time = steady_clock::now();
		last = start_time;
	}

	void stop() {
		end_time = steady_clock::now();
		total = duration<float>(end_time - start_time).count();
	}

	float split() {
		steady_clock::time_point current = steady_clock::now();
		elapsed += duration<float>(current - last).count();
		last = current;
		return elapsed;
	}
};

/* ===== Utility Functions ===== */
inline unsigned get_index(unsigned i, unsigned j, unsigned k) {
	return i + (nx * j) + (nx * ny * k);
}

inline unsigned yz_middle() {
	return 0 + (nx * int(ny / 2)) + (nx * ny * int(nz / 2));
}

inline unsigned true_middle() {
	return int(nx / 2) + int(nx * (ny / 2)) + (nx * ny * int(nz / 2));
}

#endif //REGIMPLICIT_CONSTANTS_H
