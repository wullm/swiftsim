/**
 * @triangular_cloud.cpp
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Perform triangular-shaped-cloud interpolation to determine the
 * density grid.
 */



/**
 * Function triangular_cloud_assign()
 * Perform triangular cloud interpolation
 *
 * @param rho_box Reference to array describing the density grid
 * @param k_box Reference to a vector of massive bodies
 * @param N The dimension of the grid (=width)
 * @param box_len Physical dimension of the box in Mpc
 * @param shift (default 0) Shift grid by "q" w.r.t. particles: (x,y,z) |-> (x+q,y+q,z+q)
 */

void triangular_cloud_assign(double* rho_box, std::vector<corpuscle>& bodies, int N, double box_len, double shift = 0) {
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
		int lookLftX = (int) floor((X-iX) - 1.5 + shift);
		int lookRgtX = (int) floor((X-iX) + 1.5 + shift);
		int lookLftY = (int) floor((Y-iY) - 1.5 + shift);
		int lookRgtY = (int) floor((Y-iY) + 1.5 + shift);
		int lookLftZ = (int) floor((Z-iZ) - 1.5 + shift);
		int lookRgtZ = (int) floor((Z-iZ) + 1.5 + shift);

		//Do the mass assignment
		for (int i=lookLftX; i<=lookRgtX; i++) {
			for (int j=lookLftY; j<=lookRgtY; j++) {
				for (int k=lookLftZ; k<=lookRgtZ; k++) {
					double x = abs(X - (iX+i+shift));
					double y = abs(Y - (iY+j+shift));
					double z = abs(Z - (iZ+k+shift));

					double part_x = x < 0.5 ? (0.75-x*x) : (x < 1.5 ? 0.5*pow(1.5-x, 2) : 0);
					double part_y = y < 0.5 ? (0.75-y*y) : (y < 1.5 ? 0.5*pow(1.5-y, 2) : 0);
					double part_z = z < 0.5 ? (0.75-z*z) : (z < 1.5 ? 0.5*pow(1.5-z, 2) : 0);

					rho_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] += body.mass/particle_vol * (part_x*part_y*part_z);
				}
			}
		}
	}
}
