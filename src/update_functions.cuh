#ifndef REGIMPLICIT_UPDATE_FUNCTIONS_H
#define REGIMPLICIT_UPDATE_FUNCTIONS_H

#include "constants.cuh"


/* ===== Implicit electric field updates ===== */
/* N -> N + 1/2 */
// ex = Ex + c1 * ddy * Bz - c1 * Jx
void implicit_ex_half(float* ex, const float* Ex, const float* Bz, const float* Jx) {
	for (auto k = 0; k < nz; k++) {             // [0, nz)
		for (auto j = 1; j < ny - 1; j++) {     // [1, ny - 1)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);

				auto next_y = get_index(i, j + 1, k);
				auto last_y = get_index(i, j - 1, k);

				auto temp = Ex[cur_ind] + (c1 * ddy * (Bz[next_y] - Bz[last_y])) - (c1 * Jx[cur_ind]);

				if (std::isnan(temp) || std::isinf(temp)) {
					std::cout << "NaN/Infinity in implicit_ex_half." << std::endl;
					std::exit(0);
				}

				ex[cur_ind] = temp;
			}
		}
	}
}

// ey = Ey + c1 * ddz * Bx - c1 * Jy
void implicit_ey_half(float* ey, const float* Ey, const float* Bx, const float* Jy) {
	for (auto k = 1; k < nz - 1; k++) {         // [1, nz - 1)
		for (auto j = 0; j < ny; j++) {         // [0, ny)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);

				auto next_z = get_index(i, j, k + 1);
				auto last_z = get_index(i, j, k - 1);

				auto temp = Ey[cur_ind] + (c1 * ddz * (Bx[next_z] - Bx[last_z])) - (c1 * Jy[cur_ind]);

				if (std::isnan(temp) || std::isinf(temp)) {
					std::cout << "NaN/Infinity in implicit_ey_half." << std::endl;
					std::exit(0);
				}
				ey[cur_ind] = temp;
			}
		}
	}
}

// ez = Ez + c1 * ddx * By - c1 * Jz
void implicit_ez_half(float* ez, const float* Ez, const float* By, const float* Jz) {
	for (auto k = 0; k < nz; k++) {                 // [0, nz)
		for (auto j = 0; j < ny; j++) {             // [0, ny)
			for (auto i = 1; i < nx - 1; i++) {     // [1, nx - 1)
				auto cur_ind = get_index(i, j, k);

				auto next_x = get_index(i + 1, j, k);
				auto last_x = get_index(i - 1, j, k);

				auto temp = Ez[cur_ind] + (c1 * ddx * (By[next_x] - By[last_x])) - (c1 * Jz[cur_ind]);

				if (std::isnan(temp) || std::isinf(temp)) {
					std::cout << "NaN/Infinity in implicit_ez_half." << std::endl;
					std::exit(0);
				}
				ez[cur_ind] = temp;
			}
		}
	}
}

/* N + 1/2 -> N + 1 */
// ex = Ex - c1 * ddz * By
void implicit_ex_one(float* ex, const float* Ex, const float* By) {
	for (auto k = 1; k < nz - 1; k++) {         // [1, nz - 1)
		for (auto j = 0; j < ny; j++) {         // [0, ny)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);

				auto next_z = get_index(i, j, k + 1);
				auto last_z = get_index(i, j, k - 1);

				auto temp =  Ex[cur_ind] - (c1 * ddz * (By[next_z] - By[last_z]));

				if (std::isnan(temp) || std::isinf(temp)) {
					std::cout << "NaN/Infinity in implicit_ex_one." << std::endl;
					std::exit(0);
				}

				ex[cur_ind] = temp;
			}
		}
	}
}

// ey = Ey - c1 * ddx * Bz
void implicit_ey_one(float* ey, const float* Ey, const float* Bz) {
	for (auto k = 0; k < nz; k++) {                 // [0, nz)
		for (auto j = 0; j < ny; j++) {             // [0, ny)
			for (auto i = 1; i < nx - 1; i++) {     // [1, nx - 1)
				auto cur_ind = get_index(i, j, k);

				auto next_x = get_index(i + 1, j, k);
				auto last_x = get_index(i - 1, j, k);

				auto temp = Ey[cur_ind] - (c1 * ddx * (Bz[next_x] - Bz[last_x]));

				if (std::isnan(temp) || std::isinf(temp)) {
					std::cout << "NaN/Infinity in implicit_ey_one." << std::endl;
					std::exit(0);
				}

				ey[cur_ind] = temp;
			}
		}
	}
}

// ez = Ez - c1 * ddy * Bx
void implicit_ez_one(float* ez, const float* Ez, const float* Bx) {
	for (auto k = 0; k < nz; k++) {             // [0, nz)
		for (auto j = 1; j < ny - 1; j++) {     // [1, ny - 1)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);

				auto next_y = get_index(i, j + 1, k);
				auto last_y = get_index(i, j - 1, k);

				auto temp = Ez[cur_ind] - (c1 * ddy * (Bx[next_y] - Bx[last_y]));


				if (std::isnan(temp) || std::isinf(temp)) {
					std::cout << "NaN/Infinity in implicit_ez_one." << std::endl;
					std::exit(0);
				}

				ez[cur_ind] = temp;
			}
		}
	}
}

/* ===== Explicit electric field update ===== */
/* Same for N -> N + 1/2 -> N + 1 */
void explicit_E(float* E, const float* e) {
	for (auto k = 0; k < nz; k++) {
		for (auto j = 0; j < ny; j++) {
			for (auto i = 0; i < nx; i++) {
				auto ind = get_index(i, j, k);
				auto temp = e[ind] - E[ind];

				if (std::isnan(temp) || std::isinf(temp)) {
					std::cout << "NaN/Infinity in explicit_E." << std::endl;
					std::exit(0);
				}
				E[ind] = temp;
			}
		}
	}
}

/* ===== Explicit Magnetic Field Updates ===== */
/* N -> N + 1/2 */
// Bx = Bx + c2 * ddz * ey
void explicit_bx_half(float* Bx, const float* ey) {
	for (auto k = 1; k < nz - 1; k++) {         // [1, nz - 1)
		for (auto j = 0; j < ny; j++) {         // [0, ny)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);

				auto next_z = get_index(i, j, k + 1);
				auto last_z = get_index(i, j, k - 1);

				auto temp = Bx[cur_ind] + (c2 * ddz * (ey[next_z] - ey[last_z]));

				if (std::isnan(temp) || std::isinf(temp)) {
					std::cout << "NaN/Infinity in explicit_bx_half." << std::endl;
					std::exit(0);
				}

				Bx[cur_ind] = temp;
			}
		}
	}
}

// By = By + c2 * ddx * ez
void explicit_by_half(float* By, const float* ez) {
	for (auto k = 0; k < nz; k++) {                 // [0, nz)
		for (auto j = 0; j < ny; j++) {             // [0, ny)
			for (auto i = 1; i < nx - 1; i++) {     // [1, nx - 1)
				auto cur_ind = get_index(i, j, k);

				auto next_x = get_index(i + 1, j, k);
				auto last_x = get_index(i - 1, j, k);

				auto temp = By[cur_ind] + (c2 * ddx * (ez[next_x] - ez[last_x]));

				if (std::isnan(temp) || std::isinf(temp)) {
					std::cout << "NaN/Infinity in explicit_by_half." << std::endl;
					std::exit(0);
				}

				By[cur_ind] = temp;
			}
		}
	}
}

// Bz = Bz + c2 * ddy * ex
void explicit_bz_half(float* Bz, const float* ex) {
	for (auto k = 0; k < nz; k++) {             // [0, nz)
		for (auto j = 1; j < ny - 1; j++) {     // [1, ny - 1)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);

				auto next_y = get_index(i, j + 1, k);
				auto last_y = get_index(i, j - 1, k);

				auto temp = Bz[cur_ind] + (c2 * ddy * (ex[next_y] - ex[last_y]));

				if (std::isnan(temp) || std::isinf(temp)) {
					std::cout << "NaN/Infinity in explicit_bz_half." << std::endl;
					std::exit(0);
				}

				Bz[cur_ind] = temp;
			}
		}
	}
}

/* N + 1/2 -> N + 1 */
// Bx = Bx + c2 * ddy * ez
void explicit_bx_one(float* Bx, const float* ez) {
	for (auto k = 0; k < nz; k++) {             // [0, nz)
		for (auto j = 1; j < ny - 1; j++) {     // [1, ny - 1)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);

				auto next_y = get_index(i, j + 1, k);
				auto last_y = get_index(i, j - 1, k);

				auto temp = Bx[cur_ind] - (c2 * ddy * (ez[next_y] - ez[last_y]));

				if (std::isnan(temp) || std::isinf(temp)) {
					std::cout << "NaN/Infinity in explicit_bx_one." << std::endl;
					std::exit(0);
				}

				Bx[cur_ind] = temp;
			}
		}
	}
}

// By = By - c2 * ddz * ex
void explicit_by_one(float* By, const float* ex) {
	for (auto k = 1; k < nz - 1; k++) {         // [1, nz - 1)
		for (auto j = 0; j < ny; j++) {         // [0, ny)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);

				auto next_z = get_index(i, j, k + 1);
				auto last_z = get_index(i, j, k - 1);

				auto temp = By[cur_ind] - (c2 * ddz * (ex[next_z] - ex[last_z]));

				if (std::isnan(temp) || std::isinf(temp)) {
					std::cout << "NaN/Infinity in explicit_by_one." << std::endl;
					std::exit(0);
				}

				By[cur_ind] = temp;
			}
		}
	}
}

// Bz = Bz - c2 * ddx * ey
void explicit_bz_one(float* Bz, const float* ey) {
	for (auto k = 0; k < nz; k++) {                 // [0, nz)
		for (auto j = 0; j < ny; j++) {             // [0, ny)
			for (auto i = 1; i < nx - 1; i++) {     // [1, nx - 1)
				auto cur_ind = get_index(i, j, k);

				auto next_x = get_index(i + 1, j, k);
				auto last_x = get_index(i - 1, j, k);

				auto temp = Bz[cur_ind] - (c2 * ddx * (ey[next_x] - ey[last_x]));

				if (std::isnan(temp) || std::isinf(temp)) {
					std::cout << "NaN/Infinity in explicit_bz_one." << std::endl;
					std::exit(0);
				}

				Bz[cur_ind] = temp;
			}
		}
	}
}

#endif //REGIMPLICIT_UPDATE_FUNCTIONS_H
