#ifndef REGIMPLICIT_TDMA_H
#define REGIMPLICIT_TDMA_H

#include "constants.cuh"

class Tridiagonal {
public:
	Tridiagonal() = delete;
	Tridiagonal(Tridiagonal&) = delete;

	explicit Tridiagonal(int size) :
		a(new float[size - 1]{}),
		c(new float[size - 1]{}),
		b(new float[size]{}),
		size(size)
	{}

	~Tridiagonal() {
		delete[] a;
		delete[] b;
		delete[] c;
	}

	void init(float deriv) const {
		float coeff = ((dt * dt) / (8.0f * eps0 * mu0)) * (deriv * deriv);
		for (int i = 0; i < size - 1; i++) {
			a[i] = -coeff;
			c[i] = -coeff;
		}

		for (int i = 0; i < size; i++) {
			b[i] =  0.5f + (2.0f * coeff);
		}
	}

public:
	float* a;
	float* b;
	float* c;
	int size;
};

//void onedim_solve(const Tridiagonal& t, const float* d, float* x) {
//	auto* cc = new float[nx - 1]{};
//	auto* dc = new float[nx]{};
//	/* Forward Pass */
//	// First elements of cc and dc
//	cc[0] = t.c[0] / t.b[0];
//	dc[0] = d[0] / t.b[0];
//
//	// loop until (i == nx - 2), where cc ends
//	for (int i = 1; i < nx - 1; i++) {
//		float m = 1.0f / (t.b[i] - t.a[i - 1] * cc[i - 1]);
//		cc[i] = t.c[i] * m;
//		dc[i] = (d[i] - (t.a[i - 1] * dc[i - 1])) * m;
//	}
//	// Last element of dc
//	dc[nx - 1] = (d[nx - 1] - (t.a[nx - 2] * dc[nx - 2])) / (t.b[nx - 1] - (t.a[nx - 2] * cc[nx - 2]));
//
//	/* Backwards substitution */
//	// Do last element first
//	x[nx - 1] = dc[nx - 1];
//	// Iterate backwards through current row
//	for (int i = nx - 2; i >= 0; i--) {
//		x[i] = dc[i] - (cc[i] * x[i + 1]);
//	}
//	delete[] cc;
//	delete[] dc;
//}

void ddx_solve(const Tridiagonal& t, const float* d, float* x) {
	auto* cc = new float[nx - 1]{};
	auto* dc = new float[nx]{};

	// rows of x
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			// First and last indices for current row
			auto zero_ind = get_index(     0, j, k);
			auto last_ind = get_index(nx - 1, j, k);

			/* Forward Pass */
			// First elements of cc and dc
			cc[0] = t.c[0] / t.b[0];
			dc[0] = d[zero_ind] / t.b[0];

			// loop until (i == nx - 2), where cc ends
			for (int i = 1; i < nx - 1; i++) {
				auto ind = get_index(i, j, k);

				float m = 1.0f / (t.b[i] - t.a[i - 1] * cc[i - 1]);
				cc[i] = t.c[i] * m;
				dc[i] = (d[ind] - (t.a[i - 1] * dc[i - 1])) * m;
			}
			// Last element of dc
			dc[nx - 1] = (d[last_ind] - (t.a[nx - 2] * dc[nx - 2])) / (t.b[nx - 1] - (t.a[nx - 2] * cc[nx - 2]));

			/* Backwards substitution */
			// Do last element first
			x[last_ind] = dc[nx - 1];
			// Iterate backwards through current row
			for (int i = nx - 2; i >= 0; i--) {
				auto ind = get_index(i, j, k);
				auto ind2 = get_index(i + 1, j, k);
				x[ind] = dc[i] - (cc[i] * x[ind2]);
			}
		}
	}
	delete[] cc;
	delete[] dc;
}

void ddy_solve(const Tridiagonal& t, const float* d, float* x) {
	// Create local copies once
	auto* cc = new float[t.size - 1]{};
	auto* dc = new float[t.size]{};

	// columns of y
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			// First and last indices for current row
			auto zero_ind = get_index(i,      0, k);
			auto last_ind = get_index(i, ny - 1, k);

			/* Forward Pass */
			// First elements of cc and dc
			cc[0] = t.c[0] / t.b[0];
			dc[0] = d[zero_ind] / t.b[0];

			// loop until i == n - 2, where cc ends
			for (int j = 1; j < ny - 1; j++) {
				auto ind = get_index(i, j, k);

				float m = 1.0f / (t.b[j] - (t.a[j - 1] * cc[j - 1]));
				cc[j] = t.c[j] * m;
				dc[j] = (d[ind] - (t.a[j - 1] * dc[j - 1])) * m;
			}
			// Last element of dc
			dc[ny - 1] = (d[last_ind] - (t.a[ny - 2] * dc[ny - 2])) / (t.b[ny - 1] - (t.a[ny - 2] * cc[ny - 2]));

			/* Backwards substitution */
			// Do last element first
			x[last_ind] = dc[ny - 1];
			// Iterate backwards through current row
			for (int j = ny - 2; j >= 0; j--) {
				auto ind  = get_index(i,     j, k);
				auto ind2 = get_index(i, j + 1, k);
				x[ind] = dc[j] - (cc[j] * x[ind2]);
			}
		}
	}
	delete[] cc;
	delete[] dc;
}

void ddz_solve(const Tridiagonal& t, const float* d, float* x) {
	// Create local copies once
	auto* cc = new float[nz - 1]{};
	auto* dc = new float[nz]{};

	// rows of z
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			// First and last indices for current row
			auto zero_ind = get_index(i, j, 0     );
			auto last_ind = get_index(i, j, nz - 1);

			/* Forward Pass */
			// First elements of cc and dc
			cc[0] = t.c[0] / t.b[0];
			dc[0] = d[zero_ind] / t.b[0];

			// loop until i == n - 2, where cc ends
			for (int k = 1; k < nz - 1; k++) {
				auto ind = get_index(i, j, k);

				float m = 1.0f / (t.b[k] - (t.a[k - 1] * cc[k - 1]));
				cc[k] = t.c[k] * m;
				dc[k] = (d[ind] - (t.a[k - 1] * dc[k - 1])) * m;
			}
			// Last element of dc
			dc[nz - 1] = (d[last_ind] - (t.a[nz - 2] * dc[nz - 2])) / (t.b[nz - 1] - (t.a[nz - 2] * cc[nz - 2]));

			/* Backwards substitution */
			// Do last element first
			x[last_ind] = dc[nz - 1];
			// Iterate backwards through current row
			for (int k = nz - 2; k >= 0; k--) {
				auto ind = get_index(i, j, k);
				auto ind2 = get_index(i, j, k + 1);
				x[ind] = dc[k] - (cc[k] * x[ind2]);
			}
		}
	}
	delete[] cc;
	delete[] dc;
}

#endif //REGIMPLICIT_TDMA_H
