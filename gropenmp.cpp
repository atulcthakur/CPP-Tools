#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <omp.h>  

struct Atom {
    std::string symbol;
    double x, y, z;
};

struct Vector {
    double x, y, z;
};

struct Box {
    double x, y, z;
};

struct Molecule {
    Atom O;
    Atom H1;
    Atom H2;
    Vector OH1;
    Vector OH2;
};

//
std::vector<std::vector<Atom>> read_posvelxyz(const std::string& filename, int num_frames) {
    std::ifstream file(filename);
    int num_atoms;
    double dum;
    std::string dummy_line;
    std::vector<std::vector<Atom>> all_frames;

    for (int i = 0; i < num_frames; ++i) {
        file >> num_atoms;
        std::getline(file, dummy_line);  // read the rest of the first line
        std::getline(file, dummy_line);  // read the second line into a dummy variable

        std::vector<Atom> atoms(num_atoms);
        for (Atom& atom : atoms) {
            file >> atom.symbol >> atom.x >> atom.y >> atom.z >> dum >> dum >> dum;
        }

        // Now atoms contains the atom symbols and coordinates for this frame
        //Add this frame to all_frames
        all_frames.push_back(atoms);
        }
        //Now all_frames contains all the frames from the file
        return all_frames;
    }
//

//
std::pair<std::vector<std::vector<Atom>>, std::vector<Box>> read_lammpstrj(const std::string& filename, int num_frames) {
    std::ifstream file(filename);
    std::vector<std::vector<Atom>> all_frames;
    std::vector<Box> all_boxes;

    for (int i = 0; i < num_frames; ++i) {
        std::string line;
        std::getline(file, line);  // ITEM: TIMESTEP
        std::getline(file, line);  // timestep number

        std::getline(file, line);  // ITEM: NUMBER OF ATOMS
        std::getline(file, line);
        int num_atoms = std::stoi(line);

        std::getline(file, line);  // ITEM: BOX BOUNDS pp pp pp
        Box box;
        for (int j = 0; j < 3; ++j) {
            std::getline(file, line);
            std::istringstream iss(line);
            double low, high;
            iss >> low >> high;
            if (j == 0) box.x = high - low;
            else if (j == 1) box.y = high - low;
            else if (j == 2) box.z = high - low;
        }
        all_boxes.push_back(box);

        std::getline(file, line);  // ITEM: ATOMS id type xs ys zs
        
	std::vector<Atom> atoms(num_atoms);
	for (int k = 0; k < num_atoms; ++k) {
            std::getline(file, line);
            std::istringstream iss(line);
            int id, type;
            iss >> id >> type;
            Atom& atom = atoms[id - 1];  // Use id as index (subtract 1 because indices start at 0)
            atom.x = atom.y = atom.z = 0.0;  // Initialize coordinates to zero
            iss >> atom.x >> atom.y >> atom.z;
            atom.x = atom.x*box.x;
            atom.y = atom.y*box.y;
            atom.z = atom.z*box.z;
            atom.symbol = (type == 1) ? "Na" : "Cl";  // replace with your actual atom types
        }
	all_frames.push_back(atoms);
	std::cout << "Frame" << i << "Done" << std::endl;
    }

    return {all_frames, all_boxes};
}
//

void write_frames(const std::string& filename, const std::vector<std::vector<Atom>>& all_frames) {
    int frame_index = 0;    
    std::ofstream file(filename);
    for (const auto& frame : all_frames) {
        file << frame.size() << "\n";
        file << " Atoms. Timestep: " << frame_index*10000 << "\n";  // You can replace this with your own comment
        for (const auto& atom : frame) {
            file << atom.symbol << " " << std::setprecision(8) << atom.x << " " << atom.y << " " << atom.z << "\n";
        }
        ++frame_index;
    }
}
//

Vector compute_minimum_image_vector(const Atom& atom1, const Atom& atom2, const Box& box) {
    Vector vec = {atom2.x - atom1.x, atom2.y - atom1.y, atom2.z - atom1.z};

    // Apply minimum image convention
    vec.x = vec.x - box.x * round(vec.x / box.x);
    vec.y = vec.y - box.y * round(vec.y / box.y);
    vec.z = vec.z - box.z * round(vec.z / box.z);

    return vec;
}
//

Vector compute_normalized_vector(const Atom& atom1, const Atom& atom2) {
    Vector vec = {atom2.x - atom1.x, atom2.y - atom1.y, atom2.z - atom1.z};
    double magnitude = std::sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
    return {vec.x / magnitude, vec.y / magnitude, vec.z / magnitude};
}
//

std::pair<std::vector<double>, std::vector<double>> compute_rdf(const std::vector<std::vector<Atom>>& all_frames, const std::vector<Box>& all_boxes, double max_distance, double bin_width) {
    int num_bins = static_cast<int>(max_distance / bin_width);
    std::vector<double> rdf(num_bins, 0.0);
    std::vector<double> bin_centers(num_bins, 0.0);
    int num_particles = 0;
    double density = 0.0;

    for (int frame_index = 0; frame_index < all_frames.size(); ++frame_index) {
        const auto& frame = all_frames[frame_index];
        const auto& box =  all_boxes[frame_index];
        double volume = box.x * box.y * box.z;  // Compute the volume of your system
        num_particles = frame.size();
        density = density + num_particles / volume;
    }
    double avgdensity = density / all_frames.size();
    
    #pragma omp parallel for  // Parallelize the outer for
    for (int frame_index = 0; frame_index < all_frames.size(); ++frame_index) {
        const auto& frame = all_frames[frame_index];
        const auto& box =  all_boxes[frame_index];
        //Box box =  {1.0, 1.0, 1.0};

        double volume = box.x * box.y * box.z;  // Compute the volume of your system
        num_particles = frame.size();
        // density = density + num_particles / volume;

        std::vector<double> rdf_local(num_bins, 0.0);
        for (int i = 0; i < frame.size(); ++i) {
            for (int j = i + 1; j < frame.size(); ++j) {
                Vector vec = compute_minimum_image_vector(frame[i], frame[j], box);

                double r = std::sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
                //std::cout << r << std::endl;
                if (r < max_distance) {
                    int bin = static_cast<int>(r / bin_width);
                    rdf_local[bin] += 2;
                }
            }
        }

	// Normalize the RDF and compute bin centers
	for (int i = 0; i < num_bins; ++i) {
    	    double bin_low = bin_width * i;
            double bin_high = bin_width * (i + 1);
            bin_centers[i] = (bin_low + bin_high) / 2;

            double shell_volume = 4.0 / 3.0 * M_PI * (bin_high * bin_high * bin_high - bin_low * bin_low * bin_low);
            double ideal_gas = shell_volume * avgdensity;
            rdf_local[i] /= all_frames.size() * num_particles * ideal_gas;
	
	   #pragma omp atomic
           rdf[i] += rdf_local[i];
     }
}
    return {bin_centers, rdf};
}

//
void write_rdf(const std::string& filename, const std::vector<double>& bin_centers, const std::vector<double>& rdf) {
    std::ofstream file(filename);

    // Write the bin centers and RDF values to the file
    for (size_t i = 0; i < bin_centers.size(); ++i) {
        file << std::setprecision(8) << "  " << bin_centers[i] << "   " << rdf[i] << "\n";
    }
}

//
int main(){
    std::string filename = "dump.lammpstrj";
    //std::string filename_xyz = "dump.xyz";
    int num_frames = 1000;  // replace with the actual number of frames you want to read
        
    auto [all_frames, all_boxes] = read_lammpstrj(filename, num_frames);
    //auto all_frames_xyz = read_posvelxyz(filename_xyz, num_frames);
    // Compute RDF
    double max_distance = 25;  // replace with the maximum distance of interest
    double bin_width = 0.01;  // replace with the desired bin width
    auto [bin_centers, rdf] = compute_rdf(all_frames, all_boxes, max_distance, bin_width);

    // Write RDF to file
    std::string output_filename = "rdf.dat";  // replace with your actual output filename
    write_rdf(output_filename, bin_centers, rdf);

 //   std::string output_filename = "output.xyz";  // replace with your actual output filename
 //   write_frames(output_filename, all_frames);    

    return 0;
    }
//
