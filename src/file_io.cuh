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

void process_ex(const std::string& filename, const float* Ex, const float* By) {
	std::ofstream file(filename);
	file << std::setprecision(std::numeric_limits<float>::max_digits10);

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 1; k < nz; k++) {
				auto i1 = get_index(i, j, k    );
				auto i2 = get_index(i, j, k - 1);
				file << Ex[i1] + (c1 * ddz * (By[i1] - By[i2]));
				if (i != nx - 1) {
					file << ", ";
				}
			}
			file << "\n";
		}
	}
	file.close();
}

void process_ey(const std::string& filename, const float* Ey, const float* Bz) {
	std::ofstream file(filename);
	file << std::setprecision(std::numeric_limits<float>::max_digits10);

	for (int i = 1; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				auto i1 = get_index(i    , j, k);
				auto i2 = get_index(i - 1, j, k);
				file << Ey[i1] + (c1 * ddx * (Bz[i1] - Bz[i2]));
				if (i != nx - 1) {
					file << ", ";
				}
			}
			file << "\n";
		}
	}
	file.close();
}

void process_ez(const std::string& filename, const float* Ez, const float* Bx) {
	std::ofstream file(filename);
	file << std::setprecision(std::numeric_limits<float>::max_digits10);

	for (int i = 0; i < nx; i++) {
		for (int j = 1; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				auto i1 = get_index(i, j, k);
				auto i2 = get_index(i, j - 1, k);
				file << Ez[i1] + (c1 * ddy * (Bx[i1] - Bx[i2]));
				if (i != nx - 1) {
					file << ", ";
				}
			}
			file << "\n";
		}
	}
	file.close();
}

void snapshot(int q, Efield& e, Bfield& b, Source& s) {
	std::string qs = std::to_string(q);
	process_ex("outputs/ex/t" + qs + ".csv", e.Ex, b.By);
	process_ey("outputs/ey/t" + qs + ".csv", e.Ey, b.Bz);
	process_ez("outputs/ez/t" + qs + ".csv", e.Ez, b.Bx);
}

#endif //REGIMPLICIT_FILE_IO_H
