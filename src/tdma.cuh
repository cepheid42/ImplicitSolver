#ifndef REGIMPLICIT_TDMA_H
#define REGIMPLICIT_TDMA_H

#include "constants.cuh"

class Tridiagonal {
public:
	Tridiagonal() {
		a = new float[nx - 1];
		c = new float[nx - 1];
		b = new float[nx];

		float coeff = 1.0f / (8.0f * eps0 * mu0) * (dt / dx) * (dt / dx);
		for (int i = 0; i < nx - 1; i++) {
			a[i] = -1 * coeff;
			c[i] = -1 * coeff;
			b[i] = 0.5f - coeff * (-2.0f);
		}
		b[nx - 1] = 0.5f - coeff * (-2.0f);
	}

	~Tridiagonal() {
		delete[] a;
		delete[] b;
		delete[] c;
	}

	void TDMAsolver(const float* d, float*x) const;

public:
	float* a;
	float* b;
	float* c;
};

void Tridiagonal::TDMAsolver(const float* d, float* x) const {
	// Create local copies once
	auto* cc = new float[nx - 1]();
	auto* dc = new float[nx]();

	// d and x are 3D arrays, require linear indices
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			// First and last indices for current row
			auto zero_ind = get_index(0.0f, j, k);
			auto last_ind = get_index(nx - 1, j, k);

			/* Forward Pass */
			// First elements of cc and dc
			cc[0] = c[0] / b[0];
			dc[0] = d[zero_ind] / b[0];

			// loop until i == n - 2, where cc ends
			for (int i = 0; i < nx - 1; i++) {
				auto ind = get_index(i, j, k);

				float m = 1.0f / (b[i] - a[i] * cc[i - 1]);
				cc[i] = c[i] * m;
				dc[i] = (d[ind] - a[i] * dc[i - 1]) * m;
			}
			// Last element of dc
			dc[nx - 1] = (d[last_ind] - a[nx - 1] * dc[nx - 2]) / (b[nx - 1] - a[nx - 1] * cc[nx - 2]);

			/* Backwards substitution */
			// Do last element first
			x[last_ind] = dc[nx - 1];
			// Iterate backwards through current row
			for (int i = nx - 1; i >= 0; i--) {
				auto ind = get_index(i, j, k);
				x[ind] = dc[i] - cc[i] * x[ind + 1];
			}
		}
	}
	delete[] cc;
	delete[] dc;
}

#endif //REGIMPLICIT_TDMA_H
