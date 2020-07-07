#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <cmath>
#include "constants.cuh"

using namespace std;

void TDMAsolver(const float* a, const float* b, const float* c, const float* d, float* x) {
	int n = nx;
	float cc[nx - 1] = {0.0f};
	float dc[nx]     = {0.0f};

	// Forward Pass
	cc[0] = c[0] / b[0];
	dc[0] = d[0] / b[0];

	for (int i = 1; i < n - 1; i++) {
		auto den = 1.0f / (b[i] - a[i] * cc[i - 1]);
		cc[i] = c[i] * den;
		dc[i] = (d[i] - a[i] * dc[i - 1]) * den;
	}
	dc[n - 1] = (d[n - 1] - a[n - 1] * dc[n - 2]) / (b[n - 1] - a[n - 1] * cc[n - 2]);

	// Back sub
	x[n - 1] = dc[n - 1];
	for (int i = n - 1; i >= 0; i--) {
		x[i] = dc[i] - cc[i] * x[i + 1];
	}
}

void run_loop(float* Ez, float* By, float* Jz) {
	float a[nx - 1] = {0.0f};
	float b[nx]     = {0.0f};
	float c[nx - 1] = {0.0f};
	float ez_nhalf[nx]     = {0.0f};
	float ez_nhalf_rhs[nx] = {0.0f};

	// Set arrays
	float coeff = 1.0f / (8.0f * eps0 * mu0) * (dt / dx) * (dt / dx);
	for (int i = 0; i < nx - 1; i++) {
		a[i] = -1 * coeff;
		c[i] = -1 * coeff;
		b[i] = 0.5f - coeff * (-2.0f);
	}
	b[nx - 1] = 0.5f - coeff * (-2.0f);

	// Begin time loop
	for (int q = 0; q < nt; q++) {
		// Source
		float t = (float(q) + 0.5f) * dt;
		Jz[0] = 9.89399f * sin(2.0f * pi * freq * t) * exp(-1 * pow(t - t0, n0) / (2.0f * pow(sig0, n0)));

		// Implicit update
		for (int n = 1; n < nx; n++) {
			ez_nhalf_rhs[n] = Ez[n] + (1.0f / eps0) * (0.5f * dt / dx) * (By[n] - By[n-1]) + (0.5f * dt / eps0) * Jz[n-1];
		}
		TDMAsolver(a, b, c, ez_nhalf_rhs, ez_nhalf);

		// Explicit update
		for (int n = 0; n < nx - 1; n++) {
			Ez[n] = ez_nhalf[n] - Ez[n];
			By[n] = By[n] + (1.0f / mu0) * (0.5f * dt / dx) * (ez_nhalf[n + 1] - ez_nhalf[n]);
		}

		// Implicit update (again)
		// Nothing to do here (requires more dimensions)

		if (q % 10 == 0) {
			cout << q << "/" << nt << endl;
			ofstream output("outputs/ez_q" + to_string(q) + ".csv");
			output << setprecision(numeric_limits<float>::max_digits10);

			for (int i = 0; i < nx; i++) {
				output << Ez[i] << "\n";
			}
			output.close();
		}
	}
}

int main() {
	float Jz[nx] = {0.0f};
	float Ez[nx] = {0.0f};
	float By[nx] = {0.0f};

	run_loop(Ez, By, Jz);

	return 0;
}
