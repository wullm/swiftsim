/**
 * @icgen.cpp
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Generate initial conditions for a cosmological N-body simulation.
 */

#include "icgen.h"

#include "create_grf.h"
#include "read_transfer.h"
#include "config.h"
#include "cosmo.h"
#include "H5Cpp.h"

int main() {
    //Start the clock
    auto time = std::chrono::seconds(std::time(NULL));
    long int unix_time = std::chrono::milliseconds(time).count();

    std::cout << "Welcome to the initial condition generator." << std::endl;
    std::cout << "Provide a random seed (enter 0 for current time): ";

    long int seed;
    std::cin >> seed;

    if (!seed) {
        seed = unix_time;
        std::cout << "Using timestamp as seed: " << seed << "." << std::endl << std::endl;
    } else {
        std::cout << "Using seed: " << seed << "." << std::endl << std::endl;
    }


    std::default_random_engine oracle;
    oracle.seed(seed);

    const int N = GRID_WIDTH;
    const int box_len = BOX_WIDTH; //Mpc
    const double delta_k = 2*M_PI/box_len; //Mpc^-1
    const double box_volume = box_len*box_len*box_len; //Mpc^3
    const double z_start = Z_START;

    //The primordial Gaussian random field and its Fourier transform
    double *primordial_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    fftw_complex *k_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));

    std::cout << "PHASE 0A - Cosmology check" << std::endl;
    std::cout << "1) Starting redshift z = " << z_start << ", a = " << a_scale_factor_of_z(z_start) << "." << std::endl;
    std::cout << "2) Cosmology table written to " << std::string(OUTPUT_DIR) << "cosmology.txt." << std::endl;
    std::cout << std::endl;

    test_cosmology(std::string(OUTPUT_DIR) + "cosmology.txt");

    std::cout << "PHASE 0B - Reading in transfer function files" << std::endl;

    //Neutrino density Transfer function data (loaded from CLASS)
    read_transfer(TF_ks, TF_T_rho, TF_T_theta);

	//Export transfer functions
	std::ofstream of(std::string(OUTPUT_DIR) + "transfer_functions.txt");
	of << "k(Mpc/h);T_cdm\n";

	for (int i = 0; i < TF_ks.size(); i++) {
		of << TF_ks[i] << " " << TF_T_rho[i] << std::endl;
	}

	of.close();

    std::cout << "1) Interpreted transfer functions exported to " << std::string(OUTPUT_DIR) << "transfer_functions.txt." << std::endl;
    std::cout << std::endl;



    std::cout << "PHASE 1A - Generating a Gaussian random field" << std::endl;

    //Generate a Gaussian random field
    generate_grf(oracle, k_box, N, box_len, sigma_func);

    // //Smooth it my man, then we can do the plan
    // for (int x=0; x<N; x++) {
    //     for (int y=0; y<N; y++) {
    //         for (int z=0; z<=N/2; z++) { //note that we stop at the (N/2+1)th entry
    //             double k_x = (x > N/2) ? (x - N)*delta_k : x*delta_k; //Mpc^-1
    //             double k_y = (y > N/2) ? (y - N)*delta_k : y*delta_k; //Mpc^-1
    //             double k_z = (z > N/2) ? (z - N)*delta_k : z*delta_k; //Mpc^-1
    //
    //             double k = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);
    //
    //             if (k>0) {
    //                 k_box[half_box_idx(N, x, y, z)][0] *= exp(-k*k*1*1);
    //                 k_box[half_box_idx(N, x, y, z)][1] *= exp(-k*k*1*1);
    //             }
    //         }
    //     }
    // }

    std::cout << "1) Done with generating the box (" << N << "^3 complex numbers)." << std::endl;
    std::cout << std::endl;


    //We still need to normalize the power spectrum
    const double R_filter = 8/h; //Mpc
    double integrated_sigma_8 = integrate_sigma_R(N, box_len, R_filter, sigma_func);

    std::cout << "PHASE 1B - Normalizing the random field" << std::endl;
    std::cout << "1a) Planck sigma_8 = " << sigma_8 << "." << std::endl;
    std::cout << "1b) Integrated sigma_8 = " << integrated_sigma_8 << " from unnormalized power spectrum." << std::endl;
    std::cout << "1ba) Hubble ratio: " << H_hubble_of_z(z_start)/H_hubble_of_z(0) << "." << std::endl;
    std::cout << "   See " << std::string(OUTPUT_DIR) << "sigma_8_integration.txt." << std::endl;

    //Normalize the Gaussian random field by multiplying the Fourier modes by the appropriate factor
    double global_PS_normalization = (sigma_8 / integrated_sigma_8) * (D_growth_factor(z_start) / D_growth_factor(0));

    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<=N/2; z++) {
                k_box[half_box_idx(N, x, y, z)][0] *= global_PS_normalization;
                k_box[half_box_idx(N, x, y, z)][1] *= global_PS_normalization;
            }
        }
    }

    std::cout << "2a) Growth factor at z=" << z_start << " is " << D_growth_factor(z_start) << "." << std::endl;
    std::cout << "2b) Growth factor at z=" << 0 << " is " << D_growth_factor(0) << "." << std::endl;
    std::cout << "=> The overall normalization is " << global_PS_normalization << "." << std::endl;
    std::cout << std::endl;

    std::cout << "PHASE 1C - Fourier transform of the random field" << std::endl;

    //FTT
    fftw_plan the_plan  = fftw_plan_dft_c2r_3d(N, N, N, k_box, primordial_box, FFTW_ESTIMATE);
    fftw_execute(the_plan);
    // fftw_destroy_plan(the_plan); //destroy later

    //Normalization (from Fourier conventions alone)
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<N; z++) {
                primordial_box[box_idx(N, x, y, z)] /= box_volume;
            }
        }
    }

    write_array_to_disk(std::string(OUTPUT_DIR) + "gaussian.box", primordial_box, N);

    std::cout << "1) The result has been written to " << std::string(OUTPUT_DIR) << "gaussian.box" << std::endl;
    std::cout << std::endl;


    //Next, we either load particle positions (e.g. from a glass) or
    // generate them from a grid
    bool gridgen = true;
    const long int particle_num = PARTICLE_NUM;
    std::vector<corpuscle> bodies(particle_num);
    if (gridgen) {
        std::cout << "PHASE 2A - Placing particles on a grid" << std::endl;
        for (int x=0; x<NP; x++) {
            for (int y=0; y<NP; y++) {
                for (int z=0; z<NP; z++) {
                    corpuscle body;

                    body.id = (long int) box_idx(NP, x, y, z);
                    body.X = x*(box_len/NP);
                    body.Y = y*(box_len/NP);
                    body.Z = z*(box_len/NP);

                    bodies[body.id] = body;
                }
            }
        }
        std::cout << "1) Done with placing " << NP << "^3 particles." << std::endl << std::endl;
    }

    //Compute the displacement field from the random field
    std::cout << "PHASE 2B - Compute the displacement vector field" << std::endl;
    std::cout << "1) Apply kernel to the Fourier transform of the primordial field." << std::endl;

    //Now determine the displacement field psi
    double *psi_x_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    double *psi_y_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    double *psi_z_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    fftw_complex *psi_k_x_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    fftw_complex *psi_k_y_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    fftw_complex *psi_k_z_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));

    //FTT the primordial box
    fftw_plan another_plan  = fftw_plan_dft_r2c_3d(N, N, N, primordial_box, k_box, FFTW_ESTIMATE);
    fftw_execute(another_plan);
    fftw_destroy_plan(another_plan);

    //Normalization
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<=N/2; z++) {
                k_box[half_box_idx(N, x, y, z)][0] *= box_volume/(N*N*N);
                k_box[half_box_idx(N, x, y, z)][1] *= box_volume/(N*N*N);
            }
        }
    }

    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<=N/2; z++) { //note that we stop at the (N/2+1)th entry
                double k_x = (x > N/2) ? (x - N)*delta_k : x*delta_k; //Mpc^-1
                double k_y = (y > N/2) ? (y - N)*delta_k : y*delta_k; //Mpc^-1
                double k_z = (z > N/2) ? (z - N)*delta_k : z*delta_k; //Mpc^-1

                double k = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                if (k>0) {
                    psi_k_x_box[half_box_idx(N, x, y, z)][0] = k_box[half_box_idx(N, x, y, z)][1] * k_x / (k*k);
                    psi_k_x_box[half_box_idx(N, x, y, z)][1] = -k_box[half_box_idx(N, x, y, z)][0] * k_x / (k*k);
                    psi_k_y_box[half_box_idx(N, x, y, z)][0] = k_box[half_box_idx(N, x, y, z)][1] * k_y / (k*k);
                    psi_k_y_box[half_box_idx(N, x, y, z)][1] = -k_box[half_box_idx(N, x, y, z)][0] * k_y / (k*k);
                    psi_k_z_box[half_box_idx(N, x, y, z)][0] = k_box[half_box_idx(N, x, y, z)][1] * k_z / (k*k);
                    psi_k_z_box[half_box_idx(N, x, y, z)][1] = -k_box[half_box_idx(N, x, y, z)][0] * k_z / (k*k);
                } else {
                    psi_k_x_box[half_box_idx(N, x, y, z)][0] = 0;
                    psi_k_x_box[half_box_idx(N, x, y, z)][1] = 0;
                    psi_k_y_box[half_box_idx(N, x, y, z)][0] = 0;
                    psi_k_y_box[half_box_idx(N, x, y, z)][1] = 0;
                    psi_k_z_box[half_box_idx(N, x, y, z)][0] = 0;
                    psi_k_z_box[half_box_idx(N, x, y, z)][1] = 0;
                }
            }
        }
    }

    std::cout << "2) Fourier transform back to real coordinates." << std::endl;

	//Do the IFFTs
	fftw_plan px  = fftw_plan_dft_c2r_3d(N, N, N, psi_k_x_box, psi_x_box, FFTW_ESTIMATE);
	fftw_execute(px);
	fftw_destroy_plan(px);
    fftw_free(psi_k_x_box);

	fftw_plan py  = fftw_plan_dft_c2r_3d(N, N, N, psi_k_y_box, psi_y_box, FFTW_ESTIMATE);
	fftw_execute(py);
	fftw_destroy_plan(py);
    fftw_free(psi_k_y_box);

	fftw_plan pz  = fftw_plan_dft_c2r_3d(N, N, N, psi_k_z_box, psi_z_box, FFTW_ESTIMATE);
	fftw_execute(pz);
	fftw_destroy_plan(pz);
    fftw_free(psi_k_z_box);

    //Normalization
	for (int x=0; x<N; x++) {
		for (int y=0; y<N; y++) {
			for (int z=0; z<N; z++) {
				psi_x_box[box_idx(N, x, y, z)] /= box_volume;
				psi_y_box[box_idx(N, x, y, z)] /= box_volume;
				psi_z_box[box_idx(N, x, y, z)] /= box_volume;
			}
		}
	}

    std::cout << "3) Done. The result consists of 3x" << N << "^3 real numbers." << std::endl;

    write_array_to_disk(std::string(OUTPUT_DIR) + "psi_x.box", psi_x_box, N);
	write_array_to_disk(std::string(OUTPUT_DIR) + "psi_y.box", psi_y_box, N);
	write_array_to_disk(std::string(OUTPUT_DIR) + "psi_z.box", psi_z_box, N);

    std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "psi_x.box" << std::endl;
    std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "psi_y.box" << std::endl;
    std::cout << "   Output written to " << std::string(OUTPUT_DIR) + "psi_z.box" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "PHASE 2C - Compute the velocity proportionality constant" << std::endl;

    //Compute the constant of proportionality
    double a = a_scale_factor_of_z(z_start);
    double H = H_hubble_of_z(z_start);
    double f = logarithmic_derivative_f_1(z_start);
    double dVdX = pow(a,2)*H*f;

    std::cout << "1) Expansion factor a(z) = " << a << "." << std::endl;
    std::cout << "2) Hubble rate H(z) = " << H << "." << std::endl;
    std::cout << "3) Logarithmic derivative of growth factor f(z) = " << f << "." << std::endl;
    std::cout << "4) Proportionality constant dVdX = a^2Hf = " << dVdX << "." << std::endl;
    std::cout << "  " << std::endl;

    std::cout << "PHASE 2D - Displace the particles" << std::endl;


    for (auto body : bodies) {
		double X = body.X;
		double Y = body.Y;
		double Z = body.Z;

        //Grid positions
		int iX = (int) floor(X*N/box_len);
		int iY = (int) floor(Y*N/box_len);
		int iZ = (int) floor(Z*N/box_len);

		//Intepolate the necessary fields with TSC
		float lookLength = 1.0;
		int lookLftX = (int) floor((X-iX) - lookLength);
		int lookRgtX = (int) floor((X-iX) + lookLength);
		int lookLftY = (int) floor((Y-iY) - lookLength);
		int lookRgtY = (int) floor((Y-iY) + lookLength);
		int lookLftZ = (int) floor((Z-iZ) - lookLength);
		int lookRgtZ = (int) floor((Z-iZ) + lookLength);

        //Accumulate interpolated values in psi_{xyz}
		double psi_x = 0, psi_y = 0, psi_z = 0;

		for (int i=lookLftX; i<=lookRgtX; i++) {
			for (int j=lookLftY; j<=lookRgtY; j++) {
				for (int k=lookLftZ; k<=lookRgtZ; k++) {
					//Pull the interpolated long-range force from the mesh
					double xx = abs(X - (iX+i));
					double yy = abs(Y - (iY+j));
					double zz = abs(Z - (iZ+k));

					double part_x = xx <= 1 ? 1-xx : 0;
					double part_y = yy <= 1 ? 1-yy : 0;
					double part_z = zz <= 1 ? 1-zz : 0;

					psi_x += psi_x_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
					psi_y += psi_y_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
					psi_z += psi_z_box[box_wrap_idx(N, iX+i, iY+j, iZ+k)] * (part_x*part_y*part_z);
				}
			}
		}

		body.X = X - psi_x;
		body.Y = Y - psi_y;
		body.Z = Z - psi_z;

		body.v_X = -dVdX * psi_x;
		body.v_Y = -dVdX * psi_y;
		body.v_Z = -dVdX * psi_z;

        body.mass = Omega_m*rho_crit*pow(Mpc,3)*box_volume/particle_num;

		bodies[body.id] = body;
	}

    std::cout << "1) Displaced " << particle_num << " particles." << std::endl;
    std::cout << "2) Particle mass " << bodies[0].mass << " kg." << std::endl;





    //Convert to Swift units. The below values are conversions to cgs
    double swift_unitcurrent = 1;
    double swift_unitlength = Mpc/cm;
    double swift_unitmass = 1.98848e+43;
    double swift_unittemp = 1;
    double swift_unittime = Gyr;

    //Convert the masses (lengths and times are unchanged)
    for (int i=0; i<particle_num; i++) {
        bodies[i].mass /= (swift_unitmass/kg); //from kg to swift unit
    }

    double swift_rho_crit = rho_crit * pow(Mpc,3) / (swift_unitmass/kg);
    std::cout << particle_num * bodies[0].mass / (swift_rho_crit*box_volume) << std::endl;
    std::cout << "Total mass " << particle_num * bodies[0].mass << std::endl;
    std::cout << "Volume " << box_volume << std::endl;
    std::cout << "Rho crit " << swift_rho_crit << std::endl;

    /*
        //Convert to Swift units. The below values are conversions to cgs
        double swift_unitcurrent = 1;
        double swift_unitlength = 3.08568e+24;
        double swift_unitmass = 1.98848e+43;
        double swift_unittemp = 1;
        double swift_unittime = 3.08568e+19;

        double M_hewon_to_swift = swift_unitmass/kg;
        double L_hewon_to_swift = swift_unitlength/(Mpc/cm);
        double T_hewon_to_swift = swift_unittime/Gyr;

        for (int i=0; i<PARTICLE_NUM; i++) {
            bodies[i].mass *= M_hewon_to_swift; //from kg to swift unit
            bodies[i].X *= L_hewon_to_swift; //from Mpc to swift unit
            bodies[i].Y *= L_hewon_to_swift; //from Mpc to swift unit
            bodies[i].Z *= L_hewon_to_swift; //from Mpc to swift unit
            bodies[i].v_X *= L_hewon_to_swift / T_hewon_to_swift; //from Mpc/Gyr to swift unit
            bodies[i].v_Y *= L_hewon_to_swift / T_hewon_to_swift; //from Mpc/Gyr to swift unit
            bodies[i].v_Z *= L_hewon_to_swift / T_hewon_to_swift; //from Mpc/Gyr to swift unit
        }

        double swift_rho_crit = rho_crit * M_hewon_to_swift / pow(swift_unitlength*cm,3);
        std::cout << PARTICLE_NUM * bodies[0].mass / (swift_rho_crit*box_volume*pow(L_hewon_to_swift,3)) << std::endl;
    */


    //Export the particles to an HDF5 file
    H5::H5File file(std::string(OUTPUT_DIR) + "particles.hdf5", H5F_ACC_TRUNC);

    //Write Header group
    H5::Group headerGroup(file.createGroup("/Header"));

    //Create new dataspace for attribute
    // DataSpace attr_dataspace = DataSpace(H5S_SCALAR);

// // Create new string datatype for attribute
    // H5::StrType strdatatype(H5::PredType::C_S1, 256); // of length 256 characters
    // const std::string strwritebuf ("This attribute is of type StrType");

    //Dataspace for a single number
    H5::DataSpace scalarSpace(H5S_SCALAR);
    //Data space for a row of 7 numbers (one per particle type)
    const std::size_t ptNDIMS = 1;
    hsize_t ptdims[ptNDIMS] = {7};
    H5::DataSpace row6Space(ptNDIMS, ptdims);
    //Data space for a row of 3 numbers (one per coordinate)
    const std::size_t cNDIMS = 1;
    hsize_t cdims[cNDIMS] = {3};
    H5::DataSpace row3Space(cNDIMS, cdims);

    //Write all the header attributes to the file
    H5::Attribute att0 = headerGroup.createAttribute("BoxSize", H5::PredType::NATIVE_FLOAT, row3Space);
    float boxsizes[3] = {box_len, box_len, box_len};
    att0.write(H5::PredType::NATIVE_FLOAT, boxsizes);

    H5::Attribute att1 = headerGroup.createAttribute("Dimension", H5::PredType::NATIVE_INT, scalarSpace);
    int coordinateDim = 3;
    att1.write(H5::PredType::NATIVE_INT, &coordinateDim);

    H5::Attribute att2 = headerGroup.createAttribute("Flag_Entropy_ICs", H5::PredType::NATIVE_INT, scalarSpace);
    int entropyFlag = 0;
    att2.write(H5::PredType::NATIVE_INT, &entropyFlag);

    H5::Attribute att3 = headerGroup.createAttribute("MassTable", H5::PredType::NATIVE_FLOAT, row6Space);
    float massTable[7] = {0,0,0,0,0,0,0};
    att3.write(H5::PredType::NATIVE_FLOAT, massTable);

    H5::Attribute att4 = headerGroup.createAttribute("NumFilesPerSnapshot", H5::PredType::NATIVE_INT, scalarSpace);
    int NumFilesPerSnapshot = 1;
    att4.write(H5::PredType::NATIVE_INT, &NumFilesPerSnapshot);

    H5::Attribute att5 = headerGroup.createAttribute("NumPart_Total", H5::PredType::NATIVE_INT, row6Space);
    int numparttotal[7] = {0,PARTICLE_NUM,0,0,0,0,NEUTRINO_NUM}; //number particles of type 0,1,2,3,4,5,6 (3 is dummy, doesn't exist)
    att5.write(H5::PredType::NATIVE_INT, &numparttotal);

    H5::Attribute att6 = headerGroup.createAttribute("NumPart_Total_HighWord", H5::PredType::NATIVE_INT, row6Space);
    int numparttotalhighword[7] = {0,0,0,0,0,0,0};
    att6.write(H5::PredType::NATIVE_INT, &numparttotalhighword);

    H5::Attribute att7 = headerGroup.createAttribute("Time", H5::PredType::NATIVE_FLOAT, scalarSpace);
    float Time = 0;
    att7.write(H5::PredType::NATIVE_FLOAT, &Time);




// // Set up write buffer for attribute
// // Create attribute and write to it



    //Choose particle group (dark matter is particle type 1)
    H5::Group group(file.createGroup("/PartType1"));

    float positions[3*PARTICLE_NUM];
    float internal_energy[PARTICLE_NUM];
    float masses[PARTICLE_NUM];
    long int particleIDs[PARTICLE_NUM];
    float smoothingLength[PARTICLE_NUM];
    float velocities[3*PARTICLE_NUM];

    for (int i=0; i<PARTICLE_NUM; i++) {
        positions[0+i*3] = bodies[i].X;
        positions[1+i*3] = bodies[i].Y;
        positions[2+i*3] = bodies[i].Z;

        internal_energy[i] = 1.23734;
        masses[i] = bodies[i].mass;
        particleIDs[i] = bodies[i].id;
        smoothingLength[i] = 2.4696;

        velocities[0+i*3] = bodies[i].v_X;
        velocities[1+i*3] = bodies[i].v_Y;
        velocities[2+i*3] = bodies[i].v_Z;
    }

    //std::cout << bodies[PARTICLE_NUM-1].X << " " << bodies[PARTICLE_NUM-1].Y << " " << bodies[PARTICLE_NUM-1].Z << std::endl;

    //Create two dataspaces (one 3xPARTICLE_NUM and one row 1xPARTICLE_NUM)
    const std::size_t NDIMS = 2;
    const std::size_t pNDIMS = 1;
    hsize_t dims[NDIMS] = {PARTICLE_NUM,3};
    hsize_t pdims[pNDIMS] = {PARTICLE_NUM};
    H5::DataSpace dataspace(NDIMS, dims);
    H5::DataSpace particle_dataspace(pNDIMS, pdims);

    //Write the position array to disk
    H5::DataSet dataset = group.createDataSet("Coordinates", H5::PredType::NATIVE_FLOAT, dataspace);
    dataset.write(positions, H5::PredType::NATIVE_FLOAT);

    // //Write the internal_energy array to disk
    // H5::DataSet dataset_energy = group.createDataSet("InternalEnergy", H5::PredType::NATIVE_FLOAT, particle_dataspace);
    // dataset_energy.write(internal_energy, H5::PredType::NATIVE_FLOAT);

    //Write the masses array to disk
    H5::DataSet dataset_mass = group.createDataSet("Masses", H5::PredType::NATIVE_FLOAT, particle_dataspace);
    dataset_mass.write(masses, H5::PredType::NATIVE_FLOAT);

    //Write the particleIDs array to disk
    H5::DataSet dataset_pid = group.createDataSet("ParticleIDs", H5::PredType::NATIVE_LONG, particle_dataspace);
    dataset_pid.write(particleIDs, H5::PredType::NATIVE_LONG);

    // //Write the smoothingLength array to disk
    // H5::DataSet dataset_sl = group.createDataSet("SmoothingLength", H5::PredType::NATIVE_FLOAT, particle_dataspace);
    // dataset_sl.write(smoothingLength, H5::PredType::NATIVE_FLOAT);

    //Write the velocity array to disk
    H5::DataSet dataset_v = group.createDataSet("Velocities", H5::PredType::NATIVE_FLOAT, dataspace);
    dataset_v.write(velocities, H5::PredType::NATIVE_FLOAT);


    /* Next, export the neutrinos */

    //Choose particle group (neutrinos is particle type 6)
    H5::Group group_nu(file.createGroup("/PartType6"));

    float positions_nu[3*NEUTRINO_NUM];
    float internal_energy_nu[NEUTRINO_NUM];
    float masses_nu[NEUTRINO_NUM];
    long int particleIDs_nu[NEUTRINO_NUM];
    float smoothingLength_nu[NEUTRINO_NUM];
    float velocities_nu[3*NEUTRINO_NUM];

    for (int i=0; i<NEUTRINO_NUM; i++) {
        positions_nu[0+i*3] = bodies[i].X;
        positions_nu[1+i*3] = bodies[i].Y;
        positions_nu[2+i*3] = bodies[i].Z;

        internal_energy_nu[i] = 1.23734;
        masses_nu[i] = bodies[i].mass;
        particleIDs_nu[i] = bodies[i].id + PARTICLE_NUM;
        smoothingLength_nu[i] = 2.4696;

        velocities_nu[0+i*3] = bodies[i].v_X;
        velocities_nu[1+i*3] = bodies[i].v_Y;
        velocities_nu[2+i*3] = bodies[i].v_Z;
    }

    //std::cout << bodies[PARTICLE_NUM-1].X << " " << bodies[PARTICLE_NUM-1].Y << " " << bodies[PARTICLE_NUM-1].Z << std::endl;

    //Create two dataspaces (one 3xPARTICLE_NUM and one row 1xPARTICLE_NUM)
    hsize_t dims_nu[NDIMS] = {NEUTRINO_NUM,3};
    hsize_t pdims_nu[pNDIMS] = {NEUTRINO_NUM};
    H5::DataSpace dataspace_nu(NDIMS, dims_nu);
    H5::DataSpace particle_dataspace_nu(pNDIMS, pdims_nu);

    //Write the position array to disk
    H5::DataSet dataset_nu = group_nu.createDataSet("Coordinates", H5::PredType::NATIVE_FLOAT, dataspace_nu);
    dataset_nu.write(positions_nu, H5::PredType::NATIVE_FLOAT);

    //Write the masses array to disk
    H5::DataSet dataset_mass_nu = group_nu.createDataSet("Masses", H5::PredType::NATIVE_FLOAT, particle_dataspace_nu);
    dataset_mass_nu.write(masses_nu, H5::PredType::NATIVE_FLOAT);

    //Write the particleIDs array to disk
    H5::DataSet dataset_pid_nu = group_nu.createDataSet("ParticleIDs", H5::PredType::NATIVE_LONG, particle_dataspace_nu);
    dataset_pid_nu.write(particleIDs_nu, H5::PredType::NATIVE_LONG);

    //Write the velocity array to disk
    H5::DataSet dataset_v_nu = group_nu.createDataSet("Velocities", H5::PredType::NATIVE_FLOAT, dataspace_nu);
    dataset_v_nu.write(velocities_nu, H5::PredType::NATIVE_FLOAT);






    //Write more attributes

    //Write RuntimePars group
    H5::Group runtimeGroup(file.createGroup("/RuntimePars"));

    H5::Attribute rt_att0 = runtimeGroup.createAttribute("PeriodicBoundariesOn", H5::PredType::NATIVE_INT, scalarSpace);
    int PeriodicBoundariesOn = 1;
    rt_att0.write(H5::PredType::NATIVE_INT, &PeriodicBoundariesOn);

    //Write Units group
    H5::Group unitsGroup(file.createGroup("/Units"));

    H5::Attribute units_att0 = unitsGroup.createAttribute("Unit current in cgs (U_I)", H5::PredType::NATIVE_DOUBLE, scalarSpace);
    units_att0.write(H5::PredType::NATIVE_DOUBLE, &swift_unitcurrent);

    H5::Attribute units_att1 = unitsGroup.createAttribute("Unit length in cgs (U_L)", H5::PredType::NATIVE_DOUBLE, scalarSpace);
    units_att1.write(H5::PredType::NATIVE_DOUBLE, &swift_unitlength);

    H5::Attribute units_att2 = unitsGroup.createAttribute("Unit mass in cgs (U_M)", H5::PredType::NATIVE_DOUBLE, scalarSpace);
    units_att2.write(H5::PredType::NATIVE_DOUBLE, &swift_unitmass);

    H5::Attribute units_att3 = unitsGroup.createAttribute("Unit temperature in cgs (U_T)", H5::PredType::NATIVE_DOUBLE, scalarSpace);
    units_att3.write(H5::PredType::NATIVE_DOUBLE, &swift_unittemp);

    H5::Attribute units_att4 = unitsGroup.createAttribute("Unit time in cgs (U_t)", H5::PredType::NATIVE_DOUBLE, scalarSpace);
    units_att4.write(H5::PredType::NATIVE_DOUBLE, &swift_unittime);

    file.close();



    //Export the primordial Gaussian random field to another HDF5 file
    H5::H5File grf_file(std::string(OUTPUT_DIR) + "gaussian.hdf5", H5F_ACC_TRUNC);

    //Write Header group
    H5::Group grf_headerGroup(grf_file.createGroup("/Header"));

    //Write all the header attributes to the file
    H5::Attribute grf_att0 = grf_headerGroup.createAttribute("BoxSize", H5::PredType::NATIVE_FLOAT, row3Space);
    float grf_boxsizes[3] = {box_len, box_len, box_len};
    grf_att0.write(H5::PredType::NATIVE_FLOAT, grf_boxsizes);

    H5::Attribute grf_att1 = grf_headerGroup.createAttribute("Dimension", H5::PredType::NATIVE_INT, scalarSpace);
    int grf_coordinateDim = 3;
    grf_att1.write(H5::PredType::NATIVE_INT, &grf_coordinateDim);

    //Write Field group
    H5::Group grf_fieldGroup(grf_file.createGroup("/Field"));

    //Create an N^3 dataspace
    const std::size_t grf_NDIMS = 3;
    hsize_t grf_dims[grf_NDIMS] = {N,N,N};
    H5::DataSpace grf_dataspace(grf_NDIMS, grf_dims);

    //Fourier transform back to real coordinates
    fftw_execute(the_plan);
    fftw_destroy_plan(the_plan);

    //Normalization (from Fourier conventions alone)
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<N; z++) {
                primordial_box[box_idx(N, x, y, z)] /= box_volume;
            }
        }
    }

    //Copy the Gaussian random field from FFTW format to a suitable formatted array
    float gaussian_arr[N*N*N];
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<N; z++) {
                gaussian_arr[box_idx(N, x, y, z)] = primordial_box[box_idx(N, x, y, z)];
            }
        }
    }


    //Write the gaussian field array to disk
    H5::DataSet grf_dataset = grf_fieldGroup.createDataSet("GaussianRandomField", H5::PredType::NATIVE_FLOAT, grf_dataspace);
    grf_dataset.write(gaussian_arr, H5::PredType::NATIVE_FLOAT);


    grf_file.close();
}
