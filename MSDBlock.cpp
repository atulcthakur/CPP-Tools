#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <tuple>
#include <numeric>
#include <algorithm>
#include <omp.h>

struct Vector {
    float x, y, z;
    //float x, y, z;
};

struct Atom {
    std::string symbol;
    Vector pos;
};

struct Molecule {
    Atom O;
    Atom H1;
    Atom H2;
    Vector COM;
    Vector Mu;
};

std::vector<std::vector<Atom>> read_frames(const std::string& filename, int num_atoms, int num_frames, int skip_interval) {
    std::ifstream file(filename);
    float dum;
    std::string dummy_line;
    std::vector<std::vector<Atom>> all_frames;

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    for (int i = 0; i < num_frames; ++i) {
        //file >> num_atoms;
        //std::getline(file, dummy_line);  // read the rest of the first line. Tested to work. 
        //std::getline(file, dummy_line);  // read the second line into a dummy variable
        // file.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // faster than getline. Not tested. 
        if (i % skip_interval == 0) {
            std::cout << i << "Frames Read" << std::endl;
            std::vector<Atom> atoms(num_atoms);
            for (Atom& atom : atoms) {
                file >> atom.symbol >> atom.pos.x >> atom.pos.y >> atom.pos.z >> dum >> dum >> dum;
            }

            // Now atoms contains the atom symbols and coordinates for this frame
            //Add this frame to all_frames
            all_frames.push_back(atoms);

            //if (i%5000 == 0) {
            //    std::cout << i << "Frames Read" << std::endl;
            //}
        } 
        else {
            // Skip lines for frames that are not being read
            //std::vector<Atom> atoms(num_atoms);
            for (int j = 0; j < num_atoms; ++j) {
                file >> dummy_line >> dum >> dum >> dum >> dum >> dum >> dum;
            }
        }

    }
        file.close();
        //Now all_frames contains all the frames from the file
        return all_frames;
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

void write_to_file(const std::vector<float>& array1, const std::vector<float>& array2, const std::string& filename) {
    std::ofstream file(filename);

    if (file.is_open()) {
        for (size_t i = 0; i < array1.size(); ++i) {
            file << i << " " << std::setw(12) << array1[i] << " " <<  array2[i] << "\n";
        }
        file.close();
    } else {
        std::cout << "Unable to open file: " << filename << std::endl;
    }
}

void write_to_file(const std::vector<double>& array1, const std::vector<double>& array2, const std::string& filename) {
    std::ofstream file(filename);

    if (file.is_open()) {
        for (size_t i = 0; i < array1.size(); ++i) {
            file << i << " " << std::setw(12) << array1[i] << " " <<  array2[i] << "\n";
        }
        file.close();
    } else {
        std::cout << "Unable to open file: " << filename << std::endl;
    }
}

std::tuple<std::vector<double>, std::vector<float>> calculateMSD(const std::vector<std::vector<Atom>> &trajectory){
    size_t numFrames = trajectory.size();
    size_t numParticles = trajectory[0].size();
    size_t maxLag = numFrames; // Include lag = 0

    std::vector<double> msd(maxLag, 0.0);
    std::vector<float> errorBars(maxLag, 0.0);

    // For lag = 0, MSD is always 0
    msd[0] = 0.0;
    errorBars[0] = 0.0;

    // Loop through lag times from 1 to maxLag - 1
    for (size_t lag = 1; lag < maxLag; ++lag) {
        std::vector<float> squaredDisplacements;
	//std::cout << lag << "Lag" << std::endl;

        for (size_t origin = 0; origin < numFrames - lag; ++origin) {
            for (size_t particle = 0; particle < numParticles; particle+=3) {  // Doing it only for the Oxygen atoms
                float dx = trajectory[origin + lag][particle].pos.x - trajectory[origin][particle].pos.x;
                float dy = trajectory[origin + lag][particle].pos.y - trajectory[origin][particle].pos.y;
                float dz = trajectory[origin + lag][particle].pos.z - trajectory[origin][particle].pos.z;
                float displacementSquared = dx * dx + dy * dy + dz * dz;
                squaredDisplacements.push_back(displacementSquared);
            }
        }
        
        // Compute mean MSD for this lag
        double sum = std::accumulate(squaredDisplacements.begin(), squaredDisplacements.end(), 0.0);
        msd[lag] = sum / squaredDisplacements.size();  // ntimeorigins * noxygens

        // Compute standard deviation (error bars)
        float variance = 0.0; 
        for (const auto &disp : squaredDisplacements) { 
            variance += (disp - msd[lag]) * (disp - msd[lag]); // Sum of squared deviations 
        } 
        variance /= squaredDisplacements.size(); // Normalize by total displacements
        errorBars[lag] = std::sqrt(variance);    // Square root to get standard deviation
    }

    // Return both MSD and error bars as a tuple
    return std::make_tuple(msd, errorBars);
}

// Declare a user-defined reduction for std::vector<float>
#pragma omp declare reduction(vec_float_plus : std::vector<float> : \
    std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<>())) \
    initializer(omp_priv = omp_orig)

std::tuple<std::vector<double>, std::vector<float>> calculateMSDParallel(const std::vector<std::vector<Atom>> &trajectory, size_t originsevery, size_t maxlagtime) {
    size_t numFrames = trajectory.size();
    size_t numParticles = trajectory[0].size();
    size_t maxLag = maxlagtime;

    std::vector<double> msd(maxLag, 0.0);       // To store MSD values
    std::vector<float> errorBars(maxLag, 0.0); // To store error bars

    int progress = 0; // Shared progress counter

    // For lag = 0, MSD and error bars are always 0
    msd[0] = 0.0;
    errorBars[0] = 0.0;

    // Loop through lag times in parallel
    #pragma omp parallel for schedule(dynamic, 10) shared(progress)
    for (size_t lag = 1; lag < maxLag; ++lag) {
        std::vector<float> squaredDisplacements;

        // Collect squared displacements for this lag
        // #pragma omp parallel for collapse(2) reduction(vec_float_plus : squaredDisplacements)
        for (size_t origin = 0; origin < numFrames - lag; origin += originsevery) {
            for (size_t particle = 0; particle < numParticles; particle+=3) {
                float dx = trajectory[origin + lag][particle].pos.x - trajectory[origin][particle].pos.x;
                float dy = trajectory[origin + lag][particle].pos.y - trajectory[origin][particle].pos.y;
                float dz = trajectory[origin + lag][particle].pos.z - trajectory[origin][particle].pos.z;
                float displacementSquared = dx * dx + dy * dy + dz * dz;
                squaredDisplacements.push_back(displacementSquared);
            }
        }

        // Compute mean MSD for this lag
        double sum = std::accumulate(squaredDisplacements.begin(), squaredDisplacements.end(), 0.0);
        msd[lag] = sum / squaredDisplacements.size();

        // Compute standard deviation (error bars)
        float variance = 0.0;
        //#pragma omp parallel for reduction(+:variance)
        for (size_t i = 0; i < squaredDisplacements.size(); ++i) {
            variance += (squaredDisplacements[i] - msd[lag]) * (squaredDisplacements[i] - msd[lag]);
        }
        variance /= squaredDisplacements.size();
        errorBars[lag] = std::sqrt(variance);

        // Update progress atomically
        #pragma omp atomic
        progress++;

        // Print progress every 10 completed iterations (avoid excessive printing)
        if (progress % 10 == 0) {
            #pragma omp critical
            {
                std::cout << "Progress: " << progress << "/" << maxLag << " lags processed." << std::endl;
            }
        }
    }

    // Return lag times, MSD, and error bars as a tuple
    return std::make_tuple(msd, errorBars);
}

int main() {
    std::string filename = "../posvel.xyz";
    int num_frames = 1000000;  // total number of frames
    std::vector<std::vector<Atom>> fullTrajectory = read_frames(filename, 3000, num_frames, 5);
    std::cout << "Reading Done" << std::endl;

    // Partition the trajectory into 4 equal blocks
    size_t numBlocks = 4;
    size_t blockSize = fullTrajectory.size() / 4;  // Every block should now be 50k frames

    // Vectors to store MSD for each block
    std::vector<std::vector<double>> blockMSD(numBlocks);

    omp_set_num_threads(50); 
    
    // Compute MSD for each block
    for (size_t b = 0; b < numBlocks; ++b) {
        size_t start = b * blockSize;
        size_t end = (b + 1) * blockSize;

        std::vector<std::vector<Atom>> blockTrajectory(fullTrajectory.begin() + start, fullTrajectory.begin() + end);
    
        auto [msd, error] = calculateMSDParallel(blockTrajectory, 1, 50000);
        blockMSD[b] = msd;

    }

    size_t maxLag = blockMSD[0].size(); // This should be 50k as above
    std::vector<double> averagedMSD(maxLag, 0.0);
    std::vector<double> msdStd(maxLag, 0.0);

    // For each lag time, average the MSD from all blocks and compute the standard deviation
    for (size_t lag = 0; lag < maxLag; ++lag) {
        double sum = 0.0;
        for (size_t b = 0; b < numBlocks; ++b) {
            sum += blockMSD[b][lag];
        }
        double mean = sum / numBlocks;
        averagedMSD[lag] = mean;

        double variance = 0.0;
        for (size_t b = 0; b < numBlocks; ++b) {
            variance += (blockMSD[b][lag] - mean) * (blockMSD[b][lag] - mean);
        }
        msdStd[lag] = std::sqrt(variance / numBlocks);
    }

    // Write the final averaged MSD and standard deviations to file
    write_to_file(averagedMSD, msdStd, "msd_block_avg.dat"); 

    return 0;
}
