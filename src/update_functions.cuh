#ifndef REGIMPLICIT_UPDATE_FUNCTIONS_H
#define REGIMPLICIT_UPDATE_FUNCTIONS_H

#include "e_field.cuh"
#include "b_field.cuh"
#include "sources.cuh"
#include "constants.cuh"


using float_ptr = std::unique_ptr<float[]>;

// Verified
// ex = Ex + c1 * ddy * Bz - c1 * Jx
// ey = Ey + c1 * ddz * Bx - c1 * Jy
// ez = Ez + c1 * ddx * By - c1 * Jz
void implicit_e_half(float_ptr& e, const float_ptr& E, const float_ptr& B, const float_ptr& J, float deriv) {
	for (int k = 0; k < nz; k++) {            // [0, nz)
		for (int j = 0; j < ny; j++) {        // [0, ny)
			for (int i = 1; i < nx; i++) {    // [1, nx)
				auto ind = get_index(i, j, k);
				e[ind] = E[ind] + (c1 * deriv * (B[ind] - B[ind - 1])) - (c1 * J[ind - 1]);
			}
		}
	}
}

// Verified
// ex = Ex - c1 * ddz * By
// ey = Ey - c1 * ddx * Bz
// ez = Ez - c1 * ddy * Bx
void implicit_e_one(float_ptr& e, const float_ptr& E, const float_ptr& B, float deriv) {
	for (int k = 0; k < nz; k++) {          // [0, nz)
		for (int j = 0; j < ny; j++) {      // [0, ny)
			for (int i = 1; i < nx; i++) {  // [1, nx)
				auto ind = get_index(i, j, k);
				e[ind] = E[ind] - (c1 * deriv * (B[ind] - B[ind - 1]));
			}
		}
	}
}

// Verified
// Same for N -> N + 1/2 -> N + 1
// E = e - E
void explicit_e(float_ptr& E, const float_ptr& e) {
	for (int k = 0; k < nz; k++) {          // [0, nz)
		for (int j = 0; j < ny; j++) {      // [0, ny)
			for (int i = 0; i < nx; i++) {  // [0, nx)
				auto ind = get_index(i, j, k);
				E[ind] = e[ind] - E[ind];
			}
		}
	}
}

// Verified
// Adds finite diff
// Bx = Bx + c2 * ddz * ey
// By = By + c2 * ddx * ez
// Bz = Bz + c2 * ddy * ex
void explicit_b_half(float_ptr& B, const float_ptr& e, float deriv) {
	for (int k = 0; k < nz; k++) {              // [0, nz)
		for (int j = 0; j < ny; j++) {          // [0, ny)
			for (int i = 0; i < nx - 1; i++) {  // [0, nx - 1)
				auto ind = get_index(i, j, k);

				// What happens to E[z, y, nx - 1] (the last value)?
				B[ind] = B[ind] + (c2 * deriv * (e[ind + 1] - e[ind]));
			}
		}
	}
}

// Verified
// Subtracts finite diff
// Bx = Bx - c2 * ddy * ez
// By = By - c2 * ddz * ex
// Bz = Bz - c2 * ddx * ey
void explicit_b_one(float_ptr& B, const float_ptr& e, float deriv) {
	for (int k = 0; k < nz; k++) {              // [0, nz)
		for (int j = 0; j < ny; j++) {          // [0, ny)
			for (int i = 0; i < nx - 1; i++) {  // [0, nx - 1)
				auto ind = get_index(i, j, k);

				// What happens to E[z, y, nx - 1] (the last value)?
				B[ind] = B[ind] - (c2 * deriv * (e[ind + 1] - e[ind]));
			}
		}
	}
}

#endif //REGIMPLICIT_UPDATE_FUNCTIONS_H
