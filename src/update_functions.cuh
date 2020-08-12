#ifndef REGIMPLICIT_UPDATE_FUNCTIONS_H
#define REGIMPLICIT_UPDATE_FUNCTIONS_H

#include "constants.cuh"

void implicit_ex_half(Efield& e, const Bfield& b, const Source& s) {
	for (int k = 0; k < nz; k++) {
		for (int j = 1; j < ny - 1; j++) {
			for (int i = 0; i < nx; i++) {
				auto ind = get_index(i, j, k);
				auto plus_y = get_index(i, j + 1, k);
				auto min_y = get_index(i, j - 1, k);

				auto c1 = dt / (4 * eps0 * dy);
				e.ex_rhs[ind] = e.Ex[ind] + c1 * (b.Hz[plus_y] - b.Hz[min_y]) - (dt / (2.0f * eps0)) * s.Jx[ind];
			}
		}
	}
}

void implicit_ey_half(Efield& e, const Bfield& b, const Source& s){
	for (int k = 1; k < nz - 1; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto ind = get_index(i, j, k);
				auto plus_z = get_index(i, j, k + 1);
				auto min_z = get_index(i, j, k - 1);

				auto c1 = dt / (4 * eps0 * dz);
				e.ey_rhs[ind] = e.Ey[ind] + c1 * (b.Hx[plus_z] - b.Hx[min_z]) - (dt / (2.0f * eps0)) * s.Jy[ind];
			}
		}
	}
}

void implicit_ez_half(Efield& e, const Bfield& b, const Source& s) {
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 1; i < nx - 1; i++) {
				auto ind = get_index(i, j, k);
				auto plus_x = get_index(i + 1, j, k);
				auto min_x = get_index(i - 1, j, k);

				auto c1 = dt / (4 * eps0 * dx);
				e.ez_rhs[ind] = e.Ez[ind] + c1 * (b.Hy[plus_x] - b.Hy[min_x]) - (dt / (2.0f * eps0)) * s.Jz[ind];
			}
		}
	}
}

void implicit_ex_one(Efield& e, const Bfield& b) {
	for (int k = 1; k < nz - 1; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto ind = get_index(i, j, k);
				auto plus_z = get_index(i, j, k + 1);
				auto min_z = get_index(i, j, k - 1);

				auto c1 = dt / (4 * eps0 * dz);
				e.ex_rhs[ind] = e.Ex[ind] - c1 * (b.Hy[plus_z] - b.Hy[min_z]);
			}
		}
	}
}

void implicit_ey_one(Efield& e, const Bfield& b) {
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 1; i < nx - 1; i++) {
				auto ind = get_index(i, j, k);
				auto plus_x = get_index(i + 1, j, k);
				auto min_x = get_index(i - 1, j, k);

				auto c1 = dt / (4 * eps0 * dx);
				e.ey_rhs[ind] = e.Ey[ind] - c1 * (b.Hz[plus_x] - b.Hz[min_x]);
			}
		}
	}
}

void implicit_ez_one(Efield& e, const Bfield& b) {
	for (int k = 0; k < nz; k++) {
		for (int j = 1; j < ny - 1; j++) {
			for (int i = 0; i < nx; i++) {
				auto ind = get_index(i, j, k);
				auto plus_y = get_index(i, j + 1, k);
				auto min_y = get_index(i, j - 1, k);

				auto c1 = dt / (4 * eps0 * dy);
				e.ez_rhs[ind] = e.Ez[ind] - c1 * (b.Hx[plus_y] - b.Hx[min_y]);
			}
		}
	}
}

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

void explicit_Hx_half(Bfield& b, const Efield& e) {
	for (int k = 1; k < nz - 1; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto ind = get_index(i, j, k);
				auto plus_z = get_index(i, j, k + 1);
				auto min_z = get_index(i, j, k - 1);

				auto c2 = dt / (4 * mu0 * dz);
				b.Hx[ind] = b.Hx[ind] + c2 * (e.ey[plus_z] - e.ey[min_z]); // - M;
			}
		}
	}
}

void explicit_Hy_half(Bfield& b, const Efield& e) {
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 1; i < nx - 1; i++) {
				auto ind = get_index(i, j, k);
				auto plus_x = get_index(i + 1, j, k);
				auto min_x = get_index(i - 1, j, k);

				auto c2 = dt / (4 * mu0 * dx);
				b.Hy[ind] = b.Hy[ind] + c2 * (e.ez[plus_x] - e.ez[min_x]); // - M;
			}
		}
	}
}

void explicit_Hz_half(Bfield& b, const Efield& e) {
	for (int k = 0; k < nz; k++) {
		for (int j = 1; j < ny - 1; j++) {
			for (int i = 0; i < nx; i++) {
				auto ind = get_index(i, j, k);
				auto plus_y = get_index(i, j + 1, k);
				auto min_y = get_index(i, j - 1, k);

				auto c2 = dt / (4 * mu0 * dy);
				b.Hz[ind] = b.Hz[ind] + c2 * (e.ex[plus_y] - e.ex[min_y]); // - M;
			}
		}
	}
}

void explicit_Hx_one(Bfield& b, const Efield& e){
	for (int k = 0; k < nz; k++) {
		for (int j = 1; j < ny - 1; j++) {
			for (int i = 0; i < nx; i++) {
				auto ind = get_index(i, j, k);
				auto plus_y = get_index(i, j + 1, k);
				auto min_y = get_index(i, j - 1, k);

				auto c2 = dt / (4 * mu0 * dy);
				b.Hx[ind] = b.Hx[ind] - c2 * (e.ez[plus_y] - e.ez[min_y]);
			}
		}
	}
}

void explicit_Hy_one(Bfield& b, const Efield& e) {
	for (int k = 1; k < nz - 1; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto ind = get_index(i, j, k);
				auto plus_z = get_index(i, j, k + 1);
				auto min_z = get_index(i, j, k - 1);

				auto c2 = dt / (4 * mu0 * dz);
				b.Hy[ind] = b.Hy[ind] - c2 * (e.ex[plus_z] - e.ex[min_z]);
			}
		}
	}
}

void explicit_Hz_one(Bfield& b, const Efield& e) {
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 1; i < nx - 1; i++) {
				auto ind = get_index(i, j, k);
				auto plus_x = get_index(i + 1, j, k);
				auto min_x = get_index(i - 1, j, k);

				auto c2 = dt / (4 * mu0 * dx);
				b.Hz[ind] = b.Hz[ind] - c2 * (e.ey[plus_x] - e.ey[min_x]);
			}
		}
	}
}

#endif //REGIMPLICIT_UPDATE_FUNCTIONS_H
