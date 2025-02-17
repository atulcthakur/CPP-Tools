#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <tuple>
#include <utility>

struct Vector {
    double x, y, z;
    Vector() = default;
    Vector(double xin, double yin, double zin): x(xin), y(yin), z(zin) {}

    // Substraction
    Vector operator-(const Vector &vec) const {
        return Vector(x - vec.x, y - vec.y, z - vec.z);
    }

    // Scalar multiplication
    Vector operator*(double scalar) const {
        return Vector(x * scalar, y * scalar, z * scalar);
    }

    // Scalar division
    Vector operator/(double scalar) const {
        return Vector(x / scalar, y / scalar, z / scalar);
    }

    // Product of two vectors
    Vector operator*(const Vector &vec) const {
        return Vector(x * vec.x, y * vec.y, z * vec.z);
    }

    // Addition of two vectors
    Vector operator+(const Vector &vec) const {
        return Vector(x + vec.x, y + vec.y, z + vec.z);
    }

    // Compute magnitude (norm) of vector
    double magnitude() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    // Scalar dot product
    double dot(const Vector &vec) const {
        return (x * vec.x + y * vec.y + z * vec.z);
    }
};

struct Atom {
    std::string symbol;
    Vector pos;

    Atom() = default;
    Atom(std::string symbolin, Vector posin) : symbol(symbolin), pos(posin) {}
    
    Vector operator-(const Atom &atm) const {
        return pos - atm.pos;
    }
};

struct Molecule {
    Atom O;
    Atom H1;
    Atom H2;
    //Vector OH1;
    //Vector OH2;
    Vector Mu;

    // Constructor to initialize atoms and compute displacement vectors
    Molecule() = default;
    Molecule(Atom oxygen, Atom hydrogen1, Atom hydrogen2) : O(oxygen), H1(hydrogen1), H2(hydrogen2), Mu(0, 0, 0) {}

    // Setter for Dipole Moment
    void setDipole(double q_H = 0.4238) {
        Vector OH1 = H1 - O;  // Displacement vector from O to H1
        Vector OH2 = H2 - O;  // Displacement vector from O to H2

        // Dipole moment: μ = Σ(qᵢ * rᵢ)
        Mu = (OH1 * q_H) + (OH2 * q_H);
        double magnitude = 1; // Mu.magnitude(); // Dot product
        Mu = Mu / magnitude;
    }

    Vector getCenterOfMass(double m_O = 15.999, double m_H = 1.008) const {
    // Total mass of the molecule
    double total_mass = m_O + 2 * m_H;

    // Compute center of mass
    double com_x = (O.pos.x * m_O + H1.pos.x * m_H + H2.pos.x * m_H) / total_mass;
    double com_y = (O.pos.y * m_O + H1.pos.y * m_H + H2.pos.y * m_H) / total_mass;
    double com_z = (O.pos.z * m_O + H1.pos.z * m_H + H2.pos.z * m_H) / total_mass;

    // Return as a Vector
    return Vector(com_x, com_y, com_z);
    }

};

std::vector<std::vector<Atom>> read_frames(const std::string& filename, int num_frames, int num_atoms = 192) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return {};
    }

    double dum, x, y, z;
    std::string symbol, dummy_line;
    std::vector<std::vector<Atom>> all_frames;

    for (int i = 0; i < num_frames; ++i) {
        file >> num_atoms;
        std::getline(file, dummy_line);  // read the rest of the first line
        std::getline(file, dummy_line);  // read the second line into a dummy variable
        std::vector<Atom> atoms;
        atoms.reserve(num_atoms);  // Optional, speeds up push_back()

        for (int j = 0; j < num_atoms; ++j) {
            file >> symbol >> x >> y >> z; // >> dum >> dum >> dum;
            atoms.emplace_back(symbol, Vector(x, y, z));  // Correct Atom initialization
        }

        all_frames.push_back(std::move(atoms));  // Avoid unnecessary copy

        if (i % 500 == 0) {
            std::cout << i << " Frames Read" << std::endl;
        }
    }

    return all_frames; // Nframes * natoms 
}

void write_frames(const std::string& filename, const std::vector<std::vector<Atom>>& all_frames) {
    int frame_index = 0;    
    std::ofstream file(filename);
    for (const auto& frame : all_frames) {
        file << frame.size() << "\n";
        file << " Atoms. Timestep: " << frame_index*10000 << "\n";  // You can replace this with your own comment
        for (const auto& atom : frame) {
            file << atom.symbol << " " << std::setprecision(8) << atom.pos.x << " " << atom.pos.y << " " << atom.pos.z << "\n";
        }
        ++frame_index;
    }
}


std::vector<std::vector<Molecule>> compute_molecules(const std::vector<std::vector<Atom>>& all_frames) {
    std::vector<std::vector<Molecule>> molecules_data;
    
    int nframes = all_frames.size();
    
    int natoms = all_frames[0].size();
    int num_molecules = natoms / 3;  // Assuming each molecule has 3 atoms (O, H, H)

    // Iterate over frames
    for (size_t t = 0; t < nframes; ++t) {
        std::vector<Molecule> molecules;
        molecules.reserve(num_molecules);  // Optional, preallocate memory. Speed this up.

        // Iterate over atoms (assuming each water molecule is represented by 3 atoms: O, H, H)
        for (size_t i = 0; i < natoms; i += 3) {
            // Reordered to reflect H, H, O atom order as seen in the ASE trajectory. (O, H , H)
            molecules.emplace_back(all_frames[t][i + 2],  all_frames[t][i + 1], all_frames[t][i]);  // Create Molecule directly.
            
            molecules.back().setDipole();  // Compute dipole after creating the molecule
        }

        molecules_data.emplace_back(std::move(molecules));  // Store molecules in molecules_data
    }

    return molecules_data;
}

// Function to generate k-vectors for an orthorhombic box
std::vector<Vector> generateKVectors(double Lx, double Ly, double Lz, int n_max, double k_max) {
    std::vector<Vector> kVectors;
    // Compute reciprocal lattice constants
    double bx = 2.0 * M_PI / Lx;
    double by = 2.0 * M_PI / Ly;
    double bz = 2.0 * M_PI / Lz;
    
    // Loop over integer indices
    for (int nx = -n_max; nx <= n_max; ++nx) {
        for (int ny = -n_max; ny <= n_max; ++ny) {
            for (int nz = -n_max; nz <= n_max; ++nz) {
                // Generate the k-vector
                Vector k(nx * bx, ny * by, nz * bz);
    
                // Exclude the zero vector if not needed
                if (nx == 0 && ny == 0 && nz == 0) continue;
                    
                // Only include k-vectors with magnitude less than a threshold
                double kMag = k.magnitude();
                if ( kMag > k_max) continue;
                
                kVectors.push_back(k);
            }
        }
    }
    return kVectors;
}

// Function to generate k-vectors for an orthorhombic box
std::vector<double> generateKzVectors(double Lx, double Ly, double Lz, int n_max, double k_max) {
    std::vector<double> kzVectors;

    double bz = 2.0 * M_PI / Lz;

    for (int nz = -n_max; nz <= n_max; ++nz) {
        // Generate the k-vector
        double kz;

        kz = nz * bz;
            
        if (nz == 0) continue;
                
        if (kz > k_max) continue;
                
        kzVectors.push_back(kz);
    }
    return kzVectors;
}


std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> dipole_density(const std::vector<double> &kzVectors, const std::vector<std::vector<Molecule>> &molecules_data, double binWidth, double kmax){
    //  Use kz-vectors 
    //  Project on z
    //  Compute Dipole Density
    int nframes = molecules_data.size();
    int nmols = molecules_data[0].size();

    int nBins = static_cast<int>(std::ceil(kmax / binWidth)) + 1;

    // Containers to accumulate binned results and count contributions.
    std::vector<double> dipole_density_abs(nBins, 0.0);
    std::vector<double> dipole_density_real(nBins, 0.0);
    std::vector<int> bin_counts(nBins, 0);  

    // Define the imaginary unit
    std::complex<double> I(0.0, 1.0);

    for (size_t frame = 0; frame < nframes; frame++){
        const auto &molecules = molecules_data[frame];

        std::cout << " Processing Frame:  " << frame << std::endl;

        for (const auto &k_vec : kzVectors){
            // Determine the bin index based on the magnitude of k_vec
            std::cout << " k-vector:  " << k_vec << std::endl;
            int bin_index = static_cast<int>(std::abs(k_vec) / binWidth);
            if (bin_index >= nBins) continue; // safety check, though should not occur

            for (size_t j = 0; j < nmols; j++){
                double pz_j = molecules[j].Mu.z;
                double z_j = molecules[j].O.pos.z;

                for (size_t k = j; k < nmols; k++){
                    double pz_k = molecules[k].Mu.z;
                    double z_k = molecules[k].O.pos.z;

                    // Will work without if block, but doing this halves the number of iterations over atoms. 
                    if (j == k) {
                        dipole_density_abs[bin_index] += pz_j * pz_k;
                        dipole_density_real[bin_index] += pz_j * pz_k;    
                    } else {
                        std::complex<double> phase = std::exp(-I * k_vec * (z_j - z_k)); 
                        std::complex<double> dipole_density_pair = pz_j * pz_k * phase; // j , k pair contribution 
                        // std::complex<double> dipole_density_conj_pair = std::conj(dipole_density_pair);  // k, j pair contribution
                        dipole_density_abs[bin_index] += 2.0 * std::abs(dipole_density_pair); // Magnitude simply gets doubled
                        dipole_density_real[bin_index] += 2.0 * dipole_density_pair.real(); // real part simply gets doubled
                        std::cout << " Dipole Dens: " << dipole_density_pair << "   " << j << "  " << k << std::endl;
                    }
                    bin_counts[bin_index] += (j == k) ? 1 : 2; // Wheather to add 1 or 2 depending on j == k or j != k
                }
            }
        }
        //std::cout << " Dipole Den:  " << dipole_density_abs[10] <<  "  " << dipole_density_abs[20] << "  " << dipole_density_abs[25] << std::endl;
        break;
    }

    // Normalize the binned density by the total number of contributions in each bin.
    for (size_t bin = 0; bin < nBins; bin++) {
        if (bin_counts[bin] > 0) {
            dipole_density_real[bin] = dipole_density_real[bin] / static_cast<double>(bin_counts[bin]);
            dipole_density_abs[bin] = dipole_density_abs[bin] / static_cast<double>(bin_counts[bin]);
            // std::cout << " Dipole Den:  " << dipole_density_abs[bin] << std::endl;
        }
    }

    std::vector<double> binCenters(nBins, 0.0);
    for (size_t i = 0; i < nBins; ++i) {
        binCenters[i] = (i + 0.5) * binWidth;
    }


return std::make_tuple(binCenters, dipole_density_abs, dipole_density_real);

}


void write_to_file(const std::vector<double>& corr1, const std::string& filename) {
    std::ofstream file(filename);

    if (file.is_open()) {
        for (size_t i = 0; i < corr1.size(); ++i) {
            file << i << " " << std::setprecision(14) << corr1[i] << "\n";
        }
        file.close();
    } else {
        std::cout << "Unable to open file: " << filename << std::endl;
    }
}

void write_tuple_to_file(const std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> &data, 
                        const std::string &filename){

    // Unpack the tuple
    const auto& binCenters = std::get<0>(data); 
    const auto& dipole_density_abs = std::get<1>(data);
    const auto& dipole_density_real = std::get<2>(data);

    if ( binCenters.size() != dipole_density_abs.size() || binCenters.size() != dipole_density_real.size()){
        std::cerr << "The vectors in tuple are not the same size" << std::endl;
        return;
    }

    std::ofstream file(filename);
    if (file.is_open()){
        file << "# binCenters   dipole_density_abs   dipole_density_real\n";

        for (size_t i = 0; i < binCenters.size(); i++){
            file << std::setprecision(14) << binCenters[i] << " "
            << std::setprecision(14) << dipole_density_abs[i] << " "
            << std::setprecision(14) << dipole_density_real[i] << "\n";
        }
        file.close();
    } else {
        std::cerr << "Cannot open file!" << filename << std::endl;
    }

}


int main() {
    std::string filename = "./srvmd.xyz";
    int num_frames = 1000;  // replace with the actual number of frames you want to read
    int num_atoms = 192;
    std::vector<std::vector<Atom>> all_frames = read_frames(filename, num_frames, num_atoms);

    // std::string output_filename = "output.xyz";  // replace with your actual output filename
    // write_frames(output_filename, all_frames);    
   
    std::vector<std::vector<Molecule>> molecules_data = compute_molecules(all_frames); 

    std::cout << " Computed Molecules " << std::endl;

    double Lx = 12.415238, Ly = 12.415238, Lz = 12.415238;
    int n_max = 20;
    double k_max = 10;
    std::vector<double> Kzvectors = generateKzVectors(Lx, Ly, Lz, n_max, k_max);

    std::cout << " K-Vectors Generated " << std::endl;

    double binWidth = 0.05;
    auto MuDens = dipole_density(Kzvectors, molecules_data, binWidth, k_max);

    write_tuple_to_file(MuDens, "dipoledensity.dat");
    return 0;
}