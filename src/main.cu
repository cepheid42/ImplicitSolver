#include "constants.cuh"
#include "e_field.cuh"
#include "b_field.cuh"
#include "sources.cuh"
#include "update_functions.cuh"
#include "tdma.cuh"
#include "file_io.cuh"

using namespace std;

const int step = 10;


void run_loop(Efield& e, Bfield& b, Source& s, Tridiagonal& tdm) {
	Timer update_loop_timer;
	update_loop_timer.start();

	// Begin time loop
	for (int q = 0; q < nt; q++) {
		// Source
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				auto ind = get_index(0, j, k);
				s.Jy[ind] = inc_ey(q, j);
				s.Jz[ind] = inc_ez(q, j);
			}
		}

		/* N -> N + 1/2 */
		// Implicit update
		implicit_e_half(e.ex_rhs, e.Ex, b.Bz, s.Jx, ddy); // ex = Ex + c1 * ddy * Bz - c2 * Jx
		implicit_e_half(e.ey_rhs, e.Ey, b.Bx, s.Jy, ddz); // ey = Ey + c1 * ddz * Bx - c2 * Jy
		implicit_e_half(e.ez_rhs, e.Ez, b.By, s.Jz, ddx); // ez = Ez + c1 * ddx * By - c2 * Jz

		tdm.TDMAsolver(e.ex_rhs, e.ex);
		tdm.TDMAsolver(e.ey_rhs, e.ey);
		tdm.TDMAsolver(e.ez_rhs, e.ez);

		// Explicit update
		explicit_e(e.ex, e.Ex);
		explicit_e(e.ey, e.Ey);
		explicit_e(e.ez, e.Ez);

		explicit_b(b.Bx, e.ey, ddz); // ddz, ey
		explicit_b(b.By, e.ez, ddx); // ddx, ez
		explicit_b(b.Bz, e.ex, ddy); // ddy, ex

		/* N + 1/2 -> N + 1 */
		// Implicit update
		implicit_e_one(e.ex_rhs, e.Ex, b.By, ddz); // ddz, By
		implicit_e_one(e.ey_rhs, e.Ey, b.Bz, ddx); // ddx, Bz
		implicit_e_one(e.ez_rhs, e.Ez, b.Bx, ddy); // ddy, Bx

		// Explicit update
		explicit_e(e.ex, e.Ex);
		explicit_e(e.ey, e.Ey);
		explicit_e(e.ez, e.Ez);

		explicit_b(b.Bx, e.ez, ddy); // ddy, ez
		explicit_b(b.By, e.ex, ddz); // ddz, ex
		explicit_b(b.Bz, e.ey, ddx); // ddx, ey

		if (q % step == 0) {
			cout << q << "/" << nt << ": " << update_loop_timer << "s" << endl;
			snapshot(q, e, b);
		}
	}
	update_loop_timer.stop();
	cout << "Total loop time: " << update_loop_timer.total << "s" << endl;
}

int main() {

	// todo:
	//  1.) Figure out why outputs are always 0
	//  2.) Plot slices
	//  3.) Output processing
	//  4.) Verify update functions work correctly


	save_params(step);

	Efield e;
	Bfield b;
	Source s;
	Tridiagonal tdm;

	run_loop(e, b, s, tdm);
	return 0;
}
