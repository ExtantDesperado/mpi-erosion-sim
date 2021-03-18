#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include "mpi.h"

//#include "Drop.h"
#include "ParallelHeightMap.h"
#include "ParallelDrop.h"

using namespace std::chrono;

void outputToFile(ParallelHeightMap map, std::string filename) {
    map.mapToRange(0.0, 255.0);

    std::ofstream imageFile;
    imageFile.open(filename);
    imageFile << "P2\n";
    imageFile << map.getXRes() << " " << map.getYRes() << "\n";
    imageFile << "255\n";

    for (int i = 0; i < map.getYRes(); i++) {
        for (int j = 0; j < map.getXRes(); j++) {
            imageFile << static_cast<int>(map[i][j]) << " ";
        }
        imageFile << std::endl;
    }

    imageFile.close();
}

int main(int argc, char** argv) {

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(0);

    int res = 256;
    const int NUM_DROPS = 30000;

    int numDropsPerProcess = NUM_DROPS / size;

    ParallelHeightMap map(res, res, rank);

    double d0 = MPI_Wtime();

    int count = 0;
    for (int i = 0; i < numDropsPerProcess; i++) {
        ParallelDrop d(&map, static_cast<double>(rand()) / RAND_MAX * (map.getXRes() - 1), static_cast<double>(rand()) / RAND_MAX * (map.getYRes() - 1));
        while (d.iterate());
    }

    double d1 = MPI_Wtime();

    std::cout << "Erosion time: " << d1 - d0 << std::endl;

    if (rank == 0) {
        outputToFile(map, "erosion.pgm");
    }

    MPI_Finalize();
}