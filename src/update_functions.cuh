#ifndef REGIMPLICIT_UPDATE_FUNCTIONS_H
#define REGIMPLICIT_UPDATE_FUNCTIONS_H

#include "e_field.cuh"
#include "b_field.cuh"
#include "sources.cuh"
#include "constants.cuh"


// Verified
void implicit_e_half(float* e, const float* E, const float* B, const float* J, const float deriv) {
	for (int k = 0; k < nz; k++) {            // [0, nz)
		for (int j = 0; j < ny; j++) {        // [0, ny)
			for (int i = 1; i < nx; i++) {    // [1, nx)
				auto ind = get_index(i, j, k);
				e[ind] = E[ind] + inv_eps0 * deriv * (B[ind] - B[ind - 1]) - ddeps0 * J[ind - 1];
			}
		}
	}
}

void implicit_e_one(float* e, const float* E, const float* B, const float deriv) {
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 1; i < nx; i++) {
				auto ind = get_index(i, j, k);
				e[ind] = E[ind] + inv_eps0 * deriv * (B[ind] - B[ind - 1]);
			}
		}
	}
}

// Verified
void explicit_e(const float* e, float* E) {
	for (int k = 0; k < nz; k++) {          // [0, nz)
		for (int j = 0; j < ny; j++) {      // [0, ny)
			for (int i = 0; i < nx; i++) {  // [0, nx)
				auto ind = get_index(i, j, k);
				E[ind] = e[ind] - E[ind];
			}
		}
	}
}

void explicit_b(float* B, float* e, float deriv) {
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx - 1; i++) {
				auto ind = get_index(i, j, k);
				B[ind] = B[ind] + inv_mu0 * deriv * (e[ind + 1] - e[ind]);
			}
		}
	}
}

#endif //REGIMPLICIT_UPDATE_FUNCTIONS_H
