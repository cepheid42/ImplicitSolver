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


void save_field(const std::string& filename, const float* grid) {
	std::ofstream file(filename);
	file << std::setprecision(std::numeric_limits<float>::max_digits10);

//	for (int i = 0; i < nx; i++) {
//		auto ind = i; // Just doing the [0, 0, i] row
//		file << grid[ind] << "\n";
//	}
//	file.close();


	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto ind = i + (j * ny) + (k * ny * nz);
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


void snapshot(int q, Efield& e, Bfield& b) {
	std::string qs = std::to_string(q);
	// todo: add output processing steps
	save_field("outputs/ex/t" + qs + ".csv", e.Ex);
	save_field("outputs/ey/t" + qs + ".csv", e.Ey);
	save_field("outputs/ez/t" + qs + ".csv", e.Ez);
}

#endif //REGIMPLICIT_FILE_IO_H
