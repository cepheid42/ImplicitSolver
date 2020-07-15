#ifndef REGIMPLICIT_UPDATE_FUNCTIONS_H
#define REGIMPLICIT_UPDATE_FUNCTIONS_H

#include "constants.cuh"

/* ===== Implicit electric field updates ===== */
/* N -> N + 1/2 */
// ex = Ex + c1 * ddy * Bz - c1 * Jx
void implicit_ex_half(float* ex, const float* Ex, const float* Bz, const float* Jx) {
	for (int k = 0; k < nz; k++) {
		for (int j = 1; j < ny; j++) {
			for (int i = 1; i < nx; i++) {
				auto i1 = get_index(i    , j    , k);
				auto i2 = get_index(i    , j - 1, k);
				auto i3 = get_index(i - 1, j    , k);
				float temp = Ex[i1] + (c1 * ddy * (Bz[i1] - Bz[i2])) - (c1 * Jx[i3]);
//				if (std::isnan(temp) || std::isinf(temp)) {
//					std::cout << "Nan found, not naan." << std::endl;
//				}
				ex[i1] = temp;
			}
		}
	}
}

// ey = Ey + c1 * ddz * Bx - c1 * Jy
void implicit_ey_half(float* ey, const float* Ey, const float* Bx, const float* Jy) {
	for (int k = 1; k < nz; k++) {
		for (int j = 1; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto i1 = get_index(i, j    , k    );
				auto i2 = get_index(i, j    , k - 1);
				auto i3 = get_index(i, j - 1, k    );
				float temp = Ey[i1] + (c1 * ddz * (Bx[i1] - Bx[i2])) - (c1 * Jy[i3]);
//				if (std::isnan(temp) || std::isinf(temp)) {
//					std::cout << "Nan found, not naan." << std::endl;
//				}
				ey[i1] = temp;
			}
		}
	}
}

// ez = Ez + c1 * ddx * By - c1 * Jz
void implicit_ez_half(float* ez, const float* Ez, const float* By, const float* Jz) {
	for (int k = 1; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 1; i < nx; i++) {
				auto i1 = get_index(i    , j, k    );
				auto i2 = get_index(i - 1, j, k    );
				auto i3 = get_index(i    , j, k - 1);
				float temp = Ez[i1] + (c1 * ddx * (By[i1] - By[i2])) - (c1 * Jz[i3]);
//				if (std::isnan(temp) || std::isinf(temp)) {
//					std::cout << "Nan found, not naan." << std::endl;
//				}
				ez[i1] = temp;
			}
		}
	}
}

/* N + 1/2 -> N + 1 */
// ex = Ex - c1 * ddz * By
void implicit_ex_one(float* ex, const float* Ex, const float* By) {
	for (int k = 1; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto i1 = get_index(i, j, k    );
				auto i2 = get_index(i, j, k - 1);
				float temp =  Ex[i1] - (c1 * ddz * (By[i1] - By[i2]));
//				if (std::isnan(temp) || std::isinf(temp)) {
//					std::cout << "Nan found, not naan." << std::endl;
//				}
				ex[i1] = temp;
			}
		}
	}
}

// ey = Ey - c1 * ddx * Bz
void implicit_ey_one(float* ey, const float* Ey, const float* Bz) {
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 1; i < nx; i++) {
				auto i1 = get_index(i    , j, k);
				auto i2 = get_index(i - 1, j, k);
				float temp = Ey[i1] - (c1 * ddx * (Bz[i1] - Bz[i2]));
//				if (std::isnan(temp) || std::isinf(temp)) {
//					std::cout << "Nan found, not naan." << std::endl;
//				}
				ey[i1] = temp;
			}
		}
	}
}

// ez = Ez - c1 * ddy * Bx
void implicit_ez_one(float* ez, const float* Ez, const float* Bx) {
	for (int k = 0; k < nz; k++) {
		for (int j = 1; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto i1 = get_index(i, j    , k);
				auto i2 = get_index(i, j - 1, k);
				float temp = Ez[i1] - (c1 * ddy * (Bx[i1] - Bx[i2]));
//				if (std::isnan(temp) || std::isinf(temp)) {
//					std::cout << "Nan found, not naan." << std::endl;
//				}
				ez[i1] = temp;
			}
		}
	}
}

/* ===== Explicit electric field update ===== */
/* Same for N -> N + 1/2 -> N + 1 */
void explicit_E(float* E, const float* e) {
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto ind = get_index(i, j, k);
				E[ind] = e[ind] - E[ind];
			}
		}
	}
}

/* ===== Explicit Magnetic Field Updates ===== */
/* N -> N + 1/2 */
// Bx = Bx + c2 * ddz * ey
void explicit_bx_half(float* Bx, const float* ey) {
	for (int k = 1; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto i1 = get_index(i, j, k    );
				auto i2 = get_index(i, j, k - 1);
				float temp = Bx[i1] + (c2 * ddz * (ey[i1] - ey[i2]));
//				if (std::isnan(temp) || std::isinf(temp)) {
//					std::cout << "Nan found, not naan." << std::endl;
//				}
				Bx[i1] = temp;
			}
		}
	}
}

// By = By + c2 * ddx * ez
void explicit_by_half(float* By, const float* ez) {
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 1; i < nx; i++) {
				auto i1 = get_index(i    , j, k);
				auto i2 = get_index(i - 1, j, k);
				float temp = By[i1] + (c2 * ddx * (ez[i1] - ez[i2]));
//				if (std::isnan(temp) || std::isinf(temp)) {
//					std::cout << "Nan found, not naan." << std::endl;
//				}
				By[i1] = temp;
			}
		}
	}
}

// Bz = Bz + c2 * ddy * ex
void explicit_bz_half(float* Bz, const float* ex) {
	for (int k = 0; k < nz; k++) {
		for (int j = 1; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto i1 = get_index(i, j    , k);
				auto i2 = get_index(i, j - 1, k);
				float temp = Bz[i1] + (c2 * ddy * (ex[i1] - ex[i2]));
//				if (std::isnan(temp) || std::isinf(temp)) {
//					std::cout << "Nan found, not naan." << std::endl;
//				}
				Bz[i1] = temp;
			}
		}
	}
}

/* N + 1/2 -> N + 1 */
// Bx = Bx + c2 * ddy * ez
void explicit_bx_one(float* Bx, const float* ez) {
	for (int k = 0; k < nz; k++) {
		for (int j = 1; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto i1 = get_index(i, j    , k);
				auto i2 = get_index(i, j - 1, k);
				float temp = Bx[i1] - (c2 * ddy * (ez[i1] - ez[i2]));
//				if (std::isnan(temp) || std::isinf(temp)) {
//					std::cout << "Nan found, not naan." << std::endl;
//				}
				Bx[i1] = temp;
			}
		}
	}
}

// By = By - c2 * ddz * ex
void explicit_by_one(float* By, const float* ex) {
	for (int k = 1; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto i1 = get_index(i, j, k    );
				auto i2 = get_index(i, j, k - 1);
				float temp = By[i1] - (c2 * ddz * (ex[i1] - ex[i2]));
//				if (std::isnan(temp) || std::isinf(temp)) {
//					std::cout << "Nan found, not naan." << std::endl;
//				}
				By[i1] = temp;
			}
		}
	}
}

// Bz = Bz - c2 * ddx * ey
void explicit_bz_one(float* Bz, const float* ey) {
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 1; i < nx; i++) {
				auto i1 = get_index(i    , j, k);
				auto i2 = get_index(i - 1, j, k);
				float temp = Bz[i1] - (c2 * ddx * (ey[i1] - ey[i2]));
//				if (std::isnan(temp) || std::isinf(temp)) {
//					std::cout << "Nan found, not naan." << std::endl;
//				}
				Bz[i1] = temp;
			}
		}
	}
}

#endif //REGIMPLICIT_UPDATE_FUNCTIONS_H
