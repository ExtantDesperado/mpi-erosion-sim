#pragma once
#include "PerlinNoise.h"
#include "mpi.h"
#include <iostream>
#include <algorithm>
#include <cmath>

class ParallelHeightMap {
public:

    ParallelHeightMap(int xr, int yr, int rank) {
        xRes = xr;
        yRes = yr;

        if (rank == 0) {
            MPI_Win_allocate_shared(xRes * yRes * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &map, &win);
        }
        else {
            int disp_unit;
            MPI_Aint ssize;
            MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &map, &win);
            MPI_Win_shared_query(win, 0, &ssize, &disp_unit, &map);
        }

        initHeightMap(4, 0.0, 255.0, rank);
    }

    void initHeightMap(int numOctaves, double lowerBound, double upperBound, int rank) {
        PerlinNoise noise(16);
        
        if (rank == 0) {
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
            
            for (int octave = 1; octave <= numOctaves; octave++) {

                for (int j = 0; j < yRes; j++) {
                    for (int i = 0; i < xRes; i++) {
                        double x = static_cast<double>(octave * i) / static_cast<double>(xRes) + octave;
                        double y = static_cast<double>(octave * j) / static_cast<double>(yRes) + octave;
                        map[j * xRes + i] += noise.getVal(x, y) / (1 << (octave - 1));
                    }
                }

            }

            mapToRange(0.0, 255.0);
            
            MPI_Win_unlock(0, win);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    void mapToRange(double lowerBound, double upperBound) {
        double minVal;
        double maxVal;

        for (int j = 0; j < yRes; j++) {
            for (int i = 0; i < xRes; i++) {

                if (i == 0 && j == 0) {
                    minVal = map[j * xRes + i];
                    maxVal = map[j * xRes + i];
                }
                else {
                    if (map[j * xRes + i] < minVal) {
                        minVal = map[j * xRes + i];
                    }
                    if (map[j * xRes + i] > maxVal) {
                        maxVal = map[j * xRes + i];
                    }
                }

            }
        }

        for (int j = 0; j < yRes; j++) {
            for (int i = 0; i < xRes; i++) {
                // Map to range [0, 1]
                map[j * xRes + i] = (map[j * xRes + i] - minVal) / (maxVal - minVal);
                // Map to range [lowerBound, upperBound]
                map[j * xRes + i] = map[j * xRes + i] * (upperBound - lowerBound) + lowerBound;
            }
        }
    }

    MPI_Win getWin() {
        return win;
    }

    double* operator[](int j) {
        return &map[j * xRes];
    }

    int getXRes() {
        return xRes;
    }

    int getYRes() {
        return yRes;
    }

private:
    MPI_Win win;
    double* map;
    int xRes;
    int yRes;
};