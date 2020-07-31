#ifndef REGIMPLICIT_UPDATE_FUNCTIONS_H
#define REGIMPLICIT_UPDATE_FUNCTIONS_H

#include "constants.cuh"

void punch_out(const float* grid, int i, int j, int k) {
	auto center  = get_index(i    , j    , k);
	auto plus_x  = get_index(i + 1, j    , k);
	auto minus_x = get_index(i - 1, j    , k);
	auto plus_y  = get_index(i    , j + 1, k);
	auto minus_y = get_index(i    , j - 1, k);
	auto plus_z  = get_index(i    , j    , k + 1);
	auto minus_z = get_index(i    , j    , k - 1);
	std::cout<< std::fixed
	<< "+z: " << grid[plus_z] << "\n"
	<< "+y: " << grid[plus_y] << "\n"
	<< "+x: " << grid[plus_x] << "\n"
	<< " 0: " << grid[center] << "\n"
	<< "-x: " << grid[minus_x] << "\n"
	<< "-y: " << grid[minus_y] << "\n"
	<< "-z: " << grid[minus_z] << std::endl;
}

#define check(grid, val, ind) { checkVal((grid), (val), (ind), __FILE__, __LINE__); }
void checkVal(const float* grid, float val, int ind, const char* file, int line) {
	if (val > ez0) {
		auto i = ind / (nx * ny);
		auto j = (ind / nz) % ny;
		auto k = ind % nz;

		fprintf(stderr, "%s:%d: element > ez0 at index %d (%d, %d, %d)\n", file, line, ind, i, j, k);
		punch_out(grid, i, j ,k);
		std::exit(0);
	}

	if (std::isnan(val) || std::isinf(val)) {
		auto i = ind / (nx * ny);
		auto j = (ind / nz) % ny;
		auto k = ind % nz;

		fprintf(stderr, "%s:%d: NaN/Infinity found at index %d (%d, %d, %d).\n", file, line, ind, i, j, k);
		punch_out(grid, i, j, k);
		std::exit(0);
	}
}


/* ===== Implicit electric field updates ===== */
/* N -> N + 1/2 */
// ex = Ex + c1 * dy * Bz - c1 * Jx
void implicit_ex_half(float* ex, const float* Ex, const float* Bz, const float* Jx) {
	for (auto k = 0; k < nz; k++) {             // [0, nz)
		for (auto j = 1; j < ny - 1; j++) {     // [1, ny - 1)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);
				auto next_y = get_index(i, j + 1, k);
				auto last_y = get_index(i, j - 1, k);

				auto temp = Ex[cur_ind] + (c1 * dy * (Bz[next_y] - Bz[last_y])) - (c1 * Jx[cur_ind]);

				check(Ex, temp, cur_ind)

				ex[cur_ind] = temp;
			}
		}
	}
}

// ey = Ey + c1 * dz * Bx - c1 * Jy
void implicit_ey_half(float* ey, const float* Ey, const float* Bx, const float* Jy) {
	for (auto k = 1; k < nz - 1; k++) {         // [1, nz - 1)
		for (auto j = 0; j < ny; j++) {         // [0, ny)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);

				auto next_z = get_index(i, j, k + 1);
				auto last_z = get_index(i, j, k - 1);

				auto temp = Ey[cur_ind] + (c1 * dz * (Bx[next_z] - Bx[last_z])) - (c1 * Jy[cur_ind]);

				check(Ey, temp, cur_ind)

				ey[cur_ind] = temp;
			}
		}
	}
}

// ez = Ez + c1 * dx * By - c1 * Jz
void implicit_ez_half(float* ez, const float* Ez, const float* By, const float* Jz) {
	for (auto k = 0; k < nz; k++) {                 // [0, nz)
		for (auto j = 0; j < ny; j++) {             // [0, ny)
			for (auto i = 1; i < nx - 1; i++) {     // [1, nx - 1)
				auto cur_ind = get_index(i, j, k);

				auto next_x = get_index(i + 1, j, k);
				auto last_x = get_index(i - 1, j, k);

				auto temp = Ez[cur_ind] + (c1 * dx * (By[next_x] - By[last_x])) - (c1 * Jz[cur_ind]);

				check(Ez, temp, cur_ind)

				ez[cur_ind] = temp;
			}
		}
	}
}

/* N + 1/2 -> N + 1 */
// ex = Ex - c1 * dz * By
void implicit_ex_one(float* ex, const float* Ex, const float* By) {
	for (auto k = 1; k < nz - 1; k++) {         // [1, nz - 1)
		for (auto j = 0; j < ny; j++) {         // [0, ny)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);

				auto next_z = get_index(i, j, k + 1);
				auto last_z = get_index(i, j, k - 1);

				auto temp =  Ex[cur_ind] - (c1 * dz * (By[next_z] - By[last_z]));

				check(Ex, temp, cur_ind)

				ex[cur_ind] = temp;
			}
		}
	}
}

// ey = Ey - c1 * dx * Bz
void implicit_ey_one(float* ey, const float* Ey, const float* Bz) {
	for (auto k = 0; k < nz; k++) {                 // [0, nz)
		for (auto j = 0; j < ny; j++) {             // [0, ny)
			for (auto i = 1; i < nx - 1; i++) {     // [1, nx - 1)
				auto cur_ind = get_index(i, j, k);

				auto next_x = get_index(i + 1, j, k);
				auto last_x = get_index(i - 1, j, k);

				auto temp = Ey[cur_ind] - (c1 * dx * (Bz[next_x] - Bz[last_x]));

				check(Ey, temp, cur_ind)

				ey[cur_ind] = temp;
			}
		}
	}
}

// ez = Ez - c1 * dy * Bx
void implicit_ez_one(float* ez, const float* Ez, const float* Bx) {
	for (auto k = 0; k < nz; k++) {             // [0, nz)
		for (auto j = 1; j < ny - 1; j++) {     // [1, ny - 1)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);

				auto next_y = get_index(i, j + 1, k);
				auto last_y = get_index(i, j - 1, k);

				auto temp = Ez[cur_ind] - (c1 * dy * (Bx[next_y] - Bx[last_y]));

				check(Ez, temp, cur_ind)

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

				check(E, temp, ind)

				E[ind] = temp;
			}
		}
	}
}

/* ===== Explicit Magnetic Field Updates ===== */
/* N -> N + 1/2 */
// Bx = Bx + c2 * dz * ey
void explicit_bx_half(float* Bx, const float* ey) {
	for (auto k = 1; k < nz - 1; k++) {         // [1, nz - 1)
		for (auto j = 0; j < ny; j++) {         // [0, ny)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);

				auto next_z = get_index(i, j, k + 1);
				auto last_z = get_index(i, j, k - 1);

				auto temp = Bx[cur_ind] + (c2 * dz * (ey[next_z] - ey[last_z]));

				check(Bx, temp, cur_ind)

				Bx[cur_ind] = temp;
			}
		}
	}
}

// By = By + c2 * dx * ez
void explicit_by_half(float* By, const float* ez) {
	for (auto k = 0; k < nz; k++) {                 // [0, nz)
		for (auto j = 0; j < ny; j++) {             // [0, ny)
			for (auto i = 1; i < nx - 1; i++) {     // [1, nx - 1)
				auto cur_ind = get_index(i, j, k);

				auto next_x = get_index(i + 1, j, k);
				auto last_x = get_index(i - 1, j, k);

				auto temp = By[cur_ind] + (c2 * dx * (ez[next_x] - ez[last_x]));

				check(By, temp, cur_ind)

				By[cur_ind] = temp;
			}
		}
	}
}

// Bz = Bz + c2 * dy * ex
void explicit_bz_half(float* Bz, const float* ex) {
	for (auto k = 0; k < nz; k++) {             // [0, nz)
		for (auto j = 1; j < ny - 1; j++) {     // [1, ny - 1)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);

				auto next_y = get_index(i, j + 1, k);
				auto last_y = get_index(i, j - 1, k);

				auto temp = Bz[cur_ind] + (c2 * dy * (ex[next_y] - ex[last_y]));

				check(Bz, temp, cur_ind)

				Bz[cur_ind] = temp;
			}
		}
	}
}

/* N + 1/2 -> N + 1 */
// Bx = Bx + c2 * dy * ez
void explicit_bx_one(float* Bx, const float* ez) {
	for (auto k = 0; k < nz; k++) {             // [0, nz)
		for (auto j = 1; j < ny - 1; j++) {     // [1, ny - 1)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);

				auto next_y = get_index(i, j + 1, k);
				auto last_y = get_index(i, j - 1, k);

				auto temp = Bx[cur_ind] - (c2 * dy * (ez[next_y] - ez[last_y]));

				check(Bx, temp, cur_ind)

				Bx[cur_ind] = temp;
			}
		}
	}
}

// By = By - c2 * dz * ex
void explicit_by_one(float* By, const float* ex) {
	for (auto k = 1; k < nz - 1; k++) {         // [1, nz - 1)
		for (auto j = 0; j < ny; j++) {         // [0, ny)
			for (auto i = 0; i < nx; i++) {     // [0, nx)
				auto cur_ind = get_index(i, j, k);

				auto next_z = get_index(i, j, k + 1);
				auto last_z = get_index(i, j, k - 1);

				auto temp = By[cur_ind] - (c2 * dz * (ex[next_z] - ex[last_z]));

				check(By, temp, cur_ind)

				By[cur_ind] = temp;
			}
		}
	}
}

// Bz = Bz - c2 * dx * ey
void explicit_bz_one(float* Bz, const float* ey) {
	for (auto k = 0; k < nz; k++) {                 // [0, nz)
		for (auto j = 0; j < ny; j++) {             // [0, ny)
			for (auto i = 1; i < nx - 1; i++) {     // [1, nx - 1)
				auto cur_ind = get_index(i, j, k);

				auto next_x = get_index(i + 1, j, k);
				auto last_x = get_index(i - 1, j, k);

				auto temp = Bz[cur_ind] - (c2 * dx * (ey[next_x] - ey[last_x]));

				check(Bz, temp, cur_ind)

				Bz[cur_ind] = temp;
			}
		}
	}
}

#endif //REGIMPLICIT_UPDATE_FUNCTIONS_H
