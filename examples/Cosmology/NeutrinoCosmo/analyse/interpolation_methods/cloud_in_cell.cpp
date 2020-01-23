/**
 * @cloud_in_cell.cpp
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Perform cloud-in-cell interpolation to determine the
 * density grid.
 */



/**
 * Function cloud_in_cell_assign()
 * Perform cloud-in-cell interpolation
 *
 * @param rho_box Reference to array describing the density grid
 * @param k_box Reference to a vector of massive bodies
 * @param N The dimension of the grid (=width)
 * @param box_len Physical dimension of the box in Mpc
 * @param shift (default 0) Shift grid by "q" w.r.t. particles: (x,y,z) |-> (x+q,y+q,z+q)
 */

void cloud_in_cell_assign(double* rho_box, std::vector<corpuscle>& bodies, int N, double box_len, double shift = 0) {
	for (corpuscle body : bodies) {
		double X = body.X / (box_len/N);
		double Y = body.Y / (box_len/N);
		double Z = body.Z / (box_len/N);

        double particle_vol = pow(378.9/64.0, 3);

		//Integer coordinates on the finer grid (for the long-range force)
		int iX = (int) floor(X);
		int iY = (int) floor(Y);
		int iZ = (int) floor(Z);

		//The search window with respect to the top-left-upper corner
		int lookLftX = (int) floor((X-iX) - 1 + shift);
		int lookRgtX = (int) floor((X-iX) + 1 + shift);
		int lookLftY = (int) floor((Y-iY) - 1 + shift);
		int lookRgtY = (int) floor((Y-iY) + 1 + shift);
		int lookLftZ = (int) floor((Z-iZ) - 1 + shift);
		int lookRgtZ = (int) floor((Z-iZ) + 1 + shift);

		//Do the mass assignment
		for (int i=lookLftX; i<=lookRgtX; i++) {
			for (int j=lookLftY; j<=lookRgtY; j++) {
				for (int k=lookLftZ; k<=lookRgtZ; k++) {
					double part_x = abs(X - (iX+i+shift)) <= 1 ? 1-abs(X - (iX+i+shift)) : 0;
					double part_y = abs(Y - (iY+j+shift)) <= 1 ? 1-abs(Y - (iY+j+shift)) : 0;
					double part_z = abs(Z - (iZ+k+shift)) <= 1 ? 1-abs(Z - (iZ+k+shift)) : 0;

					rho_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] += body.mass/particle_vol * (part_x*part_y*part_z);
				}
			}
		}
	}
}
