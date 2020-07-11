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


void save_field(const std::string& filename, const std::unique_ptr<float[]>& grid) {
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


void snapshot(int q, Efield& e, Bfield& b, Source& s) {
	std::string qs = std::to_string(q);
	// todo: add output processing steps
	save_field("outputs/ex/t" + qs + ".csv", e.Ex);
	save_field("outputs/ey/t" + qs + ".csv", e.Ey);
	save_field("outputs/ez/t" + qs + ".csv", e.Ez);
	save_field("outputs/jx/t" + qs + ".csv", s.Jx);
	save_field("outputs/jy/t" + qs + ".csv", s.Jy);
	save_field("outputs/jz/t" + qs + ".csv", s.Jz);
}

#endif //REGIMPLICIT_FILE_IO_H
