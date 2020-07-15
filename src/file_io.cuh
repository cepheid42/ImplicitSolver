#ifndef REGIMPLICIT_FILE_IO_H
#define REGIMPLICIT_FILE_IO_H

#include "constants.cuh"
#include "e_field.cuh"
#include "b_field.cuh"

void save_params(int step) {
	std::ofstream file("outputs/params.csv");
	file << nx << ", " << ny << ", " << nz << "\n"
	     << dx << ", " << dy << ", " << dz << "\n"
	     << nt << "\n"
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

void process_ex(const float* Ex, const float* By, const std::string& qs) {
	auto Excd = new float[nz * ny * nx]{};

	for (int k = 1; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto i1 = get_index(i, j, k    );
				auto i2 = get_index(i, j, k - 1);
				Excd[i1] = Ex[i1] + (c1 * ddz * (By[i1] - By[i2]));
			}
		}
	}

	// Save field to file
	save_field("outputs/ex/t" + qs + ".csv", Excd);
	delete[] Excd;
}

void process_ey(const float* Ey, const float* Bz, const std::string& qs) {
	auto Eycd = new float[nz * ny * nx]{};

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 1; i < nx; i++) {
				auto i1 = get_index(i    , j, k);
				auto i2 = get_index(i - 1, j, k);
				Eycd[i1] = Ey[i1] + (c1 * ddx * (Bz[i1] - Bz[i2]));

			}
		}
	}

	// save field to file
	save_field("outputs/ey/t" + qs + ".csv", Eycd);
	delete[] Eycd;
}

void process_ez(const float* Ez, const float* Bx, const std::string& qs) {
	auto Ezcd = new float[nz * ny * nx]{};

	for (int k = 0; k < nz; k++) {
		for (int j = 1; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto i1 = get_index(i, j, k);
				auto i2 = get_index(i, j - 1, k);
				Ezcd[i1] = Ez[i1] + (c1 * ddy * (Bx[i1] - Bx[i2]));
			}
		}
	}

	save_field("outputs/ez/t" + qs + ".csv", Ezcd);
	delete[] Ezcd;
}

void snapshot(int q, Efield& e, Bfield& b, Source& s) {
	std::string qs = std::to_string(q);

	process_ex(e.Ex, b.By, qs);
	process_ey(e.Ey, b.Bz, qs);
	process_ez(e.Ez, b.Bx, qs);
}

#endif //REGIMPLICIT_FILE_IO_H
