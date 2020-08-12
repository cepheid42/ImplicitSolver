#ifndef REGIMPLICIT_FILE_IO_H
#define REGIMPLICIT_FILE_IO_H

#include "constants.cuh"
#include "e_field.cuh"
#include "b_field.cuh"

void save_params(int step) {
	std::ofstream file("outputs/params.csv");
	file << std::setprecision(std::numeric_limits<float>::max_digits10);
	file << nx << ", " << ny << ", " << nz << "\n"
	     << dx << ", " << dy << ", " << dz << "\n"
	     << nt << ", " << dt << "\n"
	     << step << std::endl;
	file.close();
}

void save_field(const std::string& filename, float* grid) {
	std::ofstream file(filename);
	file << std::setprecision(std::numeric_limits<float>::max_digits10);

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto ind = get_index(i, j, k);
				file << grid[ind];
				if (i != nx - 1) {
					file << ", ";
				}
			}
			file << "\n";
		}
	}
	file.close();
}


void process_ex(const Efield& e, const Bfield& b, const std::string& qs) {
	auto Excd = new float[nz * ny * nx]{};

	for (int k = 1; k < nz - 1; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto cur_ind = get_index(i, j, k);
				auto next_z = get_index(i, j, k + 1);
				auto last_z = get_index(i, j, k - 1);

				auto c1 = dt / (4.0f * eps0 * dz);
				Excd[cur_ind] = e.Ex[cur_ind] + c1 * (b.Hy[next_z] - b.Hy[last_z]);
			}
		}
	}

	// Save field to file
	save_field("outputs/ex/t" + qs + ".csv", Excd);
	delete[] Excd;
}

void process_ey(const Efield& e, const Bfield& b, const std::string& qs) {
	auto Eycd = new float[nz * ny * nx]{};

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 1; i < nx - 1; i++) {
				auto cur_ind = get_index(i, j, k);
				auto next_x = get_index(i + 1, j, k);
				auto last_x = get_index(i - 1, j, k);

				auto c1 = dt / (4.0f * eps0 * dx);
				Eycd[cur_ind] = e.Ey[cur_ind] + c1 * (b.Hz[next_x] - b.Hz[last_x]);

			}
		}
	}

	// save field to file
	save_field("outputs/ey/t" + qs + ".csv", Eycd);
	delete[] Eycd;
}

void process_ez(const Efield& e, const Bfield& b, const std::string& qs) {
	auto Ezcd = new float[nz * ny * nx]{};

	for (int k = 0; k < nz; k++) {
		for (int j = 1; j < ny - 1; j++) {
			for (int i = 0; i < nx; i++) {
				auto cur_ind = get_index(i, j, k);
				auto next_y = get_index(i, j + 1, k);
				auto last_y = get_index(i, j - 1, k);

				auto c1 = dt / (4.0f * eps0 * dy);
				Ezcd[cur_ind] = e.Ez[cur_ind] + c1 * (b.Hx[next_y] - b.Hx[last_y]);
			}
		}
	}

	save_field("outputs/ez/t" + qs + ".csv", Ezcd);
	delete[] Ezcd;
}

void process_bx(const Efield& e, const Bfield& b, const std::string& qs) {
	auto Bxcd = new float[nz * ny * nx]{};

	for (int k = 0; k < nz; k++) {
		for (int j = 1; j < ny - 1; j++) {
			for (int i = 0; i < nx; i++) {
				auto cur_ind = get_index(i, j, k);
				auto next_y = get_index(i, j + 1, k);
				auto last_y = get_index(i, j - 1, k);

				auto c2 = dt / (4.0f * mu0 * dy);
				Bxcd[cur_ind] = b.Hx[cur_ind] + c2 * (e.Ez[next_y] - e.Ez[last_y]);
			}
		}
	}

	save_field("outputs/bx/t" + qs + ".csv", Bxcd);
	delete[] Bxcd;
}

void process_by(const Efield& e, const Bfield& b, const std::string& qs) {
	auto Bycd = new float[nz * ny * nx]{};

	for (int k = 1; k < nz - 1; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto cur_ind = get_index(i, j, k);
				auto next_z = get_index(i, j, k + 1);
				auto last_z = get_index(i, j, k - 1);

				auto c2 = dt / (4.0f * mu0 * dz);
				Bycd[cur_ind] = b.Hy[cur_ind] + c2 * (e.Ex[next_z] - e.Ex[last_z]);
			}
		}
	}

	save_field("outputs/by/t" + qs + ".csv", Bycd);
	delete[] Bycd;
}

void process_bz(const Efield& e, const Bfield& b, const std::string& qs) {
	auto Bzcd = new float[nz * ny * nx]{};

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 1; i < nx - 1; i++) {
				auto cur_ind = get_index(i, j, k);
				auto next_x = get_index(i + 1, j, k);
				auto last_x = get_index(i - 1, j, k);

				auto c2 = dt / (4.0f * mu0 * dx);
				Bzcd[cur_ind] = b.Hz[cur_ind] + c2 * (e.Ey[next_x] - e.Ey[last_x]);
			}
		}
	}

	save_field("outputs/bz/t" + qs + ".csv", Bzcd);
	delete[] Bzcd;
}

void snapshot(int q, Efield& e, Bfield& b) {
	const std::string qs = std::to_string(q);

	process_ex(e, b, qs);
	process_ey(e, b, qs);
	process_ez(e, b, qs);

	process_bx(e, b, qs);
	process_by(e, b, qs);
	process_bz(e, b, qs);
}

#endif //REGIMPLICIT_FILE_IO_H
