/**
 * @analyse.cpp
 * @author  Willem Elbers <whe@willemelbers.com>
 * @version 1.0
 * @license MIT (https://opensource.org/licenses/MIT)
 *
 * @section DESCRIPTION
 *
 * Analyse the output of cosmological N-body simulation.
 */

#include "config.h"
#include "analyse.h"
#include "field_io.h"

#include "interpolation_methods/cloud_in_cell.cpp"
#include "interpolation_methods/triangular_cloud.cpp"
#include "calc_powerspec.h"

#include "H5Cpp.h"

/* Calculate the momentum in eV, using E = a*sqrt(p^2 + m^2) ~ ap. */
double fermi_dirac_momentum(double a, double V) {
    //Some constants
    // Neutrino related constants
    const double eV = 1.602176634e-19;    // J
    const double eV_mass = 1.782662e-36;  // kg (this is eV/c^2)
    // Boltzmann's constant in units of eV/K
    const double k_b = 8.617333262145e-5;
    // The mass of the neutrino
    const double M_nu = 0.2;                // eV_mass (per species)
    const double M_nu_kg = M_nu * eV_mass;  // kg

    // The main units used in the HeWon simulation code
    const double Mpc = 3.086e22;  // m
    const double Gyr = 3.154e16;  // s
    const double km = 1000;       // m
    const double kg = 1000;       // g
    const double cm = 0.01;       // m

    // The speed of light in units of Mpc/Gyr
    const double c = 299792458 * Gyr / Mpc;
    const double cc = c*c;
    const double aa = a*a;

    double VV = V*V;

    //Calculate the length of the physical 3-velocity u=a*|dx/dt|
    double u = a*c*V/sqrt(aa*cc + aa*VV - VV);
    double gamma = 1.0/sqrt(1.0 - u*u/cc); //Lorentz factor
    double p_ph = u*gamma*M_nu_kg; //The physical 3-momentum
    double p_eV = p_ph * c/eV * pow(Mpc/Gyr,2); //in eV

    return p_eV;
}

double fermi_dirac_density(double p) {
    const double T_nu = 1.952784;  // Kelvin (today)
    // Neutrino related constants
    const double eV = 1.602176634e-19;    // J
    const double eV_mass = 1.782662e-36;  // kg (this is eV/c^2)
    // Boltzmann's constant in units of eV/K
    const double k_b = 8.617333262145e-5;

    double norm = 1.16748e+11;

    return norm * p*p / (exp(p / (k_b*T_nu)) + 1.0);
}

double fermi_dirac_density2(double p) {
    const double T_nu = 1.952784;  // Kelvin (today)
    // Neutrino related constants
    const double eV = 1.602176634e-19;    // J
    const double eV_mass = 1.782662e-36;  // kg (this is eV/c^2)
    // Boltzmann's constant in units of eV/K
    const double k_b = 8.617333262145e-5;

    double norm = 1.16748e+11;

    return norm * p*p / (exp(p / (k_b*T_nu)) + 1.0);
}

// p is in eV
double fermi_dirac_energy(double p) {
    const double M_nu = 0.2;                // eV_mass (per species)
    const double E = hypot(p, M_nu);  // eV
    return E/M_nu;
}

double fermi_dirac_ddensity(double p) {
    const double T_nu = 1.952784;  // Kelvin (today)
    // Neutrino related constants
    const double eV = 1.602176634e-19;    // J
    const double eV_mass = 1.782662e-36;  // kg (this is eV/c^2)
    // Boltzmann's constant in units of eV/K
    const double k_b = 8.617333262145e-5;
    const double kT = k_b*T_nu;

    double norm = 1.16748e+11;

    return (2*kT - exp(p/kT)*(p-2*kT))/(p*kT*(exp(p/kT)+1));
}

double fermi_dirac_exp(double p) {
    const double T_nu = 1.952784;  // Kelvin (today)
    // Neutrino related constants
    const double eV = 1.602176634e-19;    // J
    const double eV_mass = 1.782662e-36;  // kg (this is eV/c^2)
    // Boltzmann's constant in units of eV/K
    const double k_b = 8.617333262145e-5;
    const double kT = k_b*T_nu;

    double norm = 1.16748e+11;

    return exp(p/kT);
}

double fermi_dirac_kT() {
    const double T_nu = 1.952784;  // Kelvin (today)
    // Neutrino related constants
    const double eV = 1.602176634e-19;    // J
    const double eV_mass = 1.782662e-36;  // kg (this is eV/c^2)
    // Boltzmann's constant in units of eV/K
    const double k_b = 8.617333262145e-5;
    const double kT = k_b*T_nu;

    double norm = 1.16748e+11;

    return kT;
}


double sample_density(double p) {
    const double T_nu = 1.952784;  // Kelvin (today)
    // Neutrino related constants
    const double eV = 1.602176634e-19;    // J
    const double eV_mass = 1.782662e-36;  // kg (this is eV/c^2)
    // Boltzmann's constant in units of eV/K
    const double k_b = 8.617333262145e-5;

    // double norm = 8573.24;
    double norm = 1714.65;

    return norm * 1.0 / (exp(p / (5*k_b*T_nu)) + 1.0);
}

int main() {
    //Start the clock
    auto time = std::chrono::seconds(std::time(NULL));
    long int unix_time = std::chrono::milliseconds(time).count();

    bool with_selection = true;

    // std::string fname = "particles.hdf5";
    std::string fname = "box_0095_large.hdf5";
    std::string fname_ic = "particles_large.hdf5";

    std::cout << "Welcome to the N-body analyser." << std::endl;
    std::cout << "Opening " << INPUT_DIR + fname << "." << std::endl;

    std::string parttypename = "PartType6/";

    H5::H5File file(INPUT_DIR + fname, H5F_ACC_RDONLY);


    //Get the redshift from the file
    H5::Group cosmo = file.openGroup("Cosmology");
    H5::Attribute redsh_attr = cosmo.openAttribute("Redshift");
    float z;
    redsh_attr.read(H5::PredType::NATIVE_FLOAT, &z);

    std::cout << "The redshift was " << z << "." << std::endl;

    double a = 1.0/(z + 1.0);

    //Also read the starting redshift
    H5::Attribute abeg_attr = cosmo.openAttribute("a_beg");
    float a_start;
    abeg_attr.read(H5::PredType::NATIVE_FLOAT, &a_start);

    double afix = a/a_start;

    std::cout << "The starting scale factor was a_start = " << a_start << "." << std::endl;
    std::cout << "The current scale factor is a = " << a << "." << std::endl;

    // return 1;

    H5::DataSet dataset = file.openDataSet(parttypename + std::string("Coordinates"));
    H5::DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    hsize_t dims_out[2];
    int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
    int N_part = dims_out[0];
    int N_dim = dims_out[1];

    if (N_dim != 3) {
        throw std::logic_error("Expected three dimensions.");
    }

    float data_out[N_part][N_dim];

    H5::DataSpace memspace(rank, dims_out);

    dataset.read(data_out, H5::PredType::NATIVE_FLOAT, memspace, dataspace);

    std::vector<corpuscle> bodies;

    const int N = GRID_WIDTH;
    const double box_len = BOX_WIDTH; //Mpc
    const double delta_k = 2*M_PI/box_len; //Mpc^-1
    const double box_volume = box_len*box_len*box_len; //Mpc^3
    const double z_start = Z_START; //starting redshift

    for (int i=0; i<N_part; i++) {
        corpuscle body;

        body.X = data_out[i][0];
        body.Y = data_out[i][1];
        body.Z = data_out[i][2];

        bodies.push_back(body);
    }

    std::cout << "Read " << bodies.size() << " particle positions." << std::endl;

    H5::DataSet dataset2 = file.openDataSet(parttypename + std::string("Masses"));
    H5::DataSpace dataspace2 = dataset2.getSpace();
    int rank2 = dataspace2.getSimpleExtentNdims();
    hsize_t dims_out2[2];
    int ndims2 = dataspace2.getSimpleExtentDims(dims_out2, NULL);
    int N_part2 = dims_out2[0];

    if (N_part2 != N_part) {
        throw std::logic_error("Unequal number of particles.");
    }

    float data_out2[N_part2];
    dataset2.read(data_out2, H5::PredType::NATIVE_FLOAT, dataspace2, dataspace2);

    int deselected = 0;
    for (int i=0; i<N_part; i++) {
        bodies[i].mass = data_out2[i];
    }

    std::cout << deselected << " particles (" << (double) deselected/N_part*100 << "%) were deselected." << std::endl << std::endl;

    //Load the ids
    H5::DataSet dataset4 = file.openDataSet(parttypename + std::string("ParticleIDs"));
    H5::DataSpace dataspace4 = dataset4.getSpace();
    int rank4 = dataspace4.getSimpleExtentNdims();
    hsize_t dims_out4[2];
    int ndims4 = dataspace4.getSimpleExtentDims(dims_out4, NULL);
    int N_part4 = dims_out4[0];

    int data_out4[N_part4];
    dataset4.read(data_out4, H5::PredType::NATIVE_INT, dataspace4, dataspace4);

    //Find the maximum id
    int id_max = data_out4[0];
    for (int i=0; i<N_part; i++) {
        id_max = (data_out4[i] > id_max) ? data_out4[i] : id_max;
    }


    std::vector<int> particle_ids(id_max);
    for (int i=0; i<N_part; i++) {
        int id = data_out4[i];
        bodies[i].id = id;
        particle_ids[id] = i;
    }

    H5::DataSet dataset3 = file.openDataSet(parttypename + std::string("Velocities"));
    H5::DataSpace dataspace3 = dataset3.getSpace();
    int rank3 = dataspace3.getSimpleExtentNdims();
    hsize_t dims_out3[2];
    int ndims3 = dataspace3.getSimpleExtentDims(dims_out3, NULL);
    int N_part3 = dims_out3[0];
    int N_dims3 = dims_out3[1];

    if (N_part3 != N_part) {
        throw std::logic_error("Unequal number of particles.");
    }

    dataset3.read(data_out, H5::PredType::NATIVE_FLOAT, dataspace3, dataspace3);

    for (int i=0; i<N_part; i++) {
        double v_X = data_out[i][0];
        double v_Y = data_out[i][1];
        double v_Z = data_out[i][2];

        bodies[i].v_X = v_X;
        bodies[i].v_Y = v_Y;
        bodies[i].v_Z = v_Z;
        bodies[i].V = sqrt(v_X*v_X + v_Y*v_Y + v_Z*v_Z);
    }

    //Determine the slowest and fastest speed
    double Vmin = bodies[0].V;
    double Vmax = bodies[0].V;
    double Vavg = 0;
    double Msum = 0;

    for (int i=1; i<N_part; i++) {
        Vmin = (bodies[i].V < Vmin) ? bodies[i].V : Vmin;
        Vmax = (bodies[i].V > Vmax) ? bodies[i].V : Vmax;
        Vavg += bodies[i].V*bodies[i].mass;
        Msum += bodies[i].mass;
    }



    Vavg /= Msum;

    Vmin=0.580054;
    Vmax=124.693;

    std::cout << std::endl;
    std::cout << "V_min = " << Vmin << std::endl;
    std::cout << "V_max = " << Vmax << std::endl;
    std::cout << "V_avg = " << Vavg << std::endl;
    std::cout << std::endl;

    //Load the ic velocities
    H5::H5File file_ic(INPUT_DIR + fname_ic, H5F_ACC_RDONLY);
    H5::DataSet dataset_ic = file_ic.openDataSet(parttypename + std::string("Velocities"));

    H5::DataSpace dataspace_ic = dataset_ic.getSpace();
    int rank_ic = dataspace_ic.getSimpleExtentNdims();
    hsize_t dims_out_ic[2];
    int ndims_ic = dataspace_ic.getSimpleExtentDims(dims_out_ic, NULL);
    int N_part_ic = dims_out_ic[0];
    int N_dim_ic = dims_out_ic[1];

    if (N_dim_ic != 3) {
        throw std::logic_error("Expected three dimensions.");
    } else if (N_part_ic != N_part) {
        throw std::logic_error("Unequal numbers of particles.");
    }

    // float data_out_ic[N_part][N_dim];

    std::cout << N_part_ic << std::endl;
    std::cout << N_dim_ic << std::endl;

    //Reuse data_out
    dataset_ic.read(data_out, H5::PredType::NATIVE_FLOAT, dataspace_ic, dataspace_ic);

    for (int i=0; i<N_part; i++) {
        int id = particle_ids[64*64*64+i];
        // int id = particle_ids[i];

        double v_X = data_out[i][0];
        double v_Y = data_out[i][1];
        double v_Z = data_out[i][2];

        bodies[id].ic_v_X = v_X;
        bodies[id].ic_v_Y = v_Y;
        bodies[id].ic_v_Z = v_Z;
        bodies[id].ic_V = sqrt(v_X*v_X + v_Y*v_Y + v_Z*v_Z);
    }


    deselected = 0;

    for (int i=0; i<N_part; i++) {
        corpuscle body = bodies[i];
        double p_eV = fermi_dirac_momentum(a, a*body.V);
        double p_eV_ic = fermi_dirac_momentum(a_start, a_start*body.ic_V);
        double f = fermi_dirac_density(p_eV);
        double f_ic = fermi_dirac_density(p_eV_ic);
        double g_ic = sample_density(p_eV_ic);

        double beta = 0.1/hypot(0.1, a*body.V);
    }

    std::cout << deselected << " particles (" << (double) deselected/N_part*100 << "%) were deselected." << std::endl << std::endl;

    int annoying = 0;
    for (int i=0; i<N_part; i++) {
        corpuscle body = bodies[i];
        double p_eV = fermi_dirac_momentum(a, a*body.V);
        double p_eV_ic = fermi_dirac_momentum(a_start, a_start*body.ic_V);
        double f = fermi_dirac_density(p_eV);
        double f_ic = fermi_dirac_density(p_eV_ic);
        double g_ic = sample_density(p_eV_ic);
    }

    std::cout << (double) annoying/N_part << " annoying particles." << std::endl;

    const int nbins = 100;
    int count[nbins];
    double sep[nbins];

    for (int i=0; i<nbins; i++) {
        count[i] = 0;
        sep[i] = i*(Vmax-Vmin)/nbins;
    }

    for (int i=0; i<N_part; i++) {
        double X = bodies[i].X / 376. * 64.;
        double Y = bodies[i].Y / 376. * 64.;
        double Z = bodies[i].Z / 376. * 64.;

        // if (bodies[i].mass == 0) {
            int bin = (int) floor(nbins*(bodies[i].ic_V - Vmin)/(Vmax-Vmin));
            count[bin]++;
        // }
    }

    for (int i=0; i<nbins; i++) {
        std::cout << i << " " << sep[i] << " " << count[i] << std::endl;
    }

    std::cout << std::endl;

    double *rho_box = (double*) fftw_malloc(sizeof(double)*N*N*N);
    fftw_complex *rho_k_box = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2+1));
	fftw_plan p_rho  = fftw_plan_dft_r2c_3d(N, N, N, rho_box, rho_k_box, FFTW_ESTIMATE);


    //Initialize the box
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<N; z++) {
                rho_box[box_idx(N, x, y, z)] = 0;
            }
        }
    }

    triangular_cloud_assign(rho_box, bodies, N, box_len);
    // cloud_in_cell_assign(rho_box, bodies, N, box_len);

    write_array_to_disk(OUTPUT_DIR + std::string("rho.box"), rho_box, N);
    std::cout << "Exported density estimate to " << OUTPUT_DIR + std::string("rho.box") << "." << std::endl;


    double *marker_box = (double*) fftw_malloc(sizeof(double)*N*N*N);

    //Initialize the box
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<N; z++) {
                marker_box[box_idx(N, x, y, z)] = 0;
            }
        }
    }

    for (int i=0; i<N_part; i++) {
        bodies[i].mass = 1.0;
    }

    triangular_cloud_assign(marker_box, bodies, N, box_len);


    write_array_to_disk(OUTPUT_DIR + std::string("marker.box"), marker_box, N);
    std::cout << "Exported density estimate to " << OUTPUT_DIR + std::string("marker.box") << "." << std::endl;



    double *two_box = (double*) fftw_malloc(sizeof(double)*N*N*N);

    for (int i=0; i<N_part; i++) {
        corpuscle body = bodies[i];
        double p_eV = fermi_dirac_momentum(a, a*body.V);
        double p_eV_ic = fermi_dirac_momentum(a_start, a_start*body.ic_V);
        double f = fermi_dirac_density(p_eV);
        double f_ic = fermi_dirac_density(p_eV_ic);
        double g_ic = sample_density(p_eV_ic);

        bodies[i].mass = pow(p_eV - p_eV_ic, 2);
    }

    triangular_cloud_assign(two_box, bodies, N, box_len);

    write_array_to_disk(OUTPUT_DIR + std::string("two.box"), two_box, N);
    std::cout << "Exported density estimate to " << OUTPUT_DIR + std::string("two.box") << "." << std::endl;



    //Compute the ratio
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<N; z++) {
                marker_box[box_idx(N, x, y, z)] = (1-marker_box[box_idx(N, x, y, z)])*rho_box[box_idx(N, x, y, z)]/two_box[box_idx(N, x, y, z)];
            }
        }
    }

    write_array_to_disk(OUTPUT_DIR + std::string("ratio.box"), marker_box, N);
    std::cout << "Exported density estimate to " << OUTPUT_DIR + std::string("ratio.box") << "." << std::endl;



    //Calculate the power spectrum
	int bins = 40;
	double *k_in_bins = (double*) malloc(sizeof(double)*bins);
	double *power_in_bins = (double*) malloc(sizeof(double)*bins);
	int *obs_in_bins = (int*) malloc(sizeof(int)*bins);

	for (int i=0; i<bins; i++) {
		k_in_bins[i] = 0;
		power_in_bins[i] = 0;
		obs_in_bins[i] = 0;
	}

    fftw_execute(p_rho);

    //Normalization
    for (int x=0; x<N; x++) {
        for (int y=0; y<N; y++) {
            for (int z=0; z<=N/2; z++) {
                rho_k_box[half_box_idx(N, x, y, z)][0] *= box_volume/(N*N*N);
                rho_k_box[half_box_idx(N, x, y, z)][1] *= box_volume/(N*N*N);
            }
        }
    }


	calc_powerspec(N, box_len, rho_k_box, bins, k_in_bins, power_in_bins, obs_in_bins, CLOUD_IN_CELL);

    for (int i=0; i<bins; i++) {
        if (obs_in_bins[i] > 0) {
            std::cout << i << " " << k_in_bins[i] << " " << power_in_bins[i] << std::endl;
        }
	}

    for (int i=0; i<N_part; i++) {


        if (bodies[i].id==262147) {
            double vx = bodies[i].v_X;
            double vy = bodies[i].v_Y;
            double vz = bodies[i].v_Z;
            double V = a*bodies[i].V;
            double p = fermi_dirac_momentum(a, V);
            double E = fermi_dirac_energy(p);
            std::cout << "...." << a*vx << " " << a*vy << " " << a*vz << std::endl;
            std::cout << "...." << a*bodies[i].V << " " << (E-0.2)/0.2 << std::endl;
        }
    }

    std::cout << "a = " << a << std::endl;

    //
    std::cout << std::endl << std::endl;
    int counter = 0;
    for (auto body : bodies) {
        double p_eV = fermi_dirac_momentum(a, a*body.V);
        double p_eV_ic = fermi_dirac_momentum(a_start, a_start*body.ic_V);
        double f = fermi_dirac_density(p_eV);
        double f_ic = fermi_dirac_density(p_eV_ic);
        double g_ic = sample_density(p_eV_ic);

    }
    std::cout << "Count: " << counter << std::endl;
    std::cout << std::endl;


    //Determine the most and least massive particles
    double Mmin = bodies[0].mass;
    double Mmax = bodies[0].mass;
    double I = 0;

    for (int i=1; i<N_part; i++) {
        Mmin = (bodies[i].mass < Mmin) ? bodies[i].mass : Mmin;
        Mmax = (bodies[i].mass > Mmax) ? bodies[i].mass : Mmax;
        I += 0.5*pow(bodies[i].mass,2)/(N_part*37.180455);
    }


    std::cout << "M_min = " << Mmin << std::endl;
    std::cout << "M_max = " << Mmax << std::endl;
    std::cout << "M_sum = " << Msum << std::endl;
    std::cout << "I = " << I << std::endl;

    double a1 = 0.0243902;
    double V = 0.63251794636;
    double p = fermi_dirac_momentum(a1,V);
    double E = fermi_dirac_energy(p);

    printf("%.10e \t%.10e\n", p, E);
}
