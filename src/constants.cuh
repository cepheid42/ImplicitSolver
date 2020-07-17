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

//const float freq = 14.0e6f;
//const float lambda = c0 / freq;
//const float T    = 1.0f / freq;
//const float t0   = 3.0f * T;
//const float sig0 = 1.0f * T;
//const float n0   = 2.0f;
//const float ez0 = 140.0f/ 1.0f; // V/m

/* ====== Spatial Constants ====== */
const int x_resolution = 16;
const int num_wavelengths = 6;

const int nx = 50;
const int ny = 50;
const int nz = 50;

const float dx = 0.02; // 2 mm
const float dy = dx;
const float dz = dx;

/* ====== Temporal Constants ====== */
// Courant number
const float cfl = 1.0f;
const float dt = dx / (c0 * 1.732050807569f);  // ~= dx / (c * sqrt(3))

/*
 * The CFL number (S_c) is defined as (dt / dt_cfl) in the paper
 * Since S_c == 1.0, dt == dt_cfl
 *
 * if S_c == 4.0 (later in paper), then what would dt be? Is dt_cfl a constant number?
 */

const int nt = 38944;  // 150 ns / dt

const float tau = 1.5e-10f; // 150 ps
const float t0 = 3 * tau;

/* ====== Derivative Coefficients ====== */
const float ddx = 1.0f / dx;
const float ddy = 1.0f / dy;
const float ddz = 1.0f / dz;

const float c1 = dt / (2.0f * eps0);
const float c2 = dt / (2.0f * mu0);

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

#endif //REGIMPLICIT_CONSTANTS_H
