#pragma once
#include "PerlinNoise.h"

class HeightMap {
public:

    HeightMap(int xr, int yr) {
        xRes = xr;
        yRes = yr;
        initHeightMap(4, 0.0, 255.0);
    }

    void initHeightMap(int numOctaves, double lowerBound, double upperBound) {
        PerlinNoise noise(16);

        map = new double* [yRes];
        for (int j = 0; j < yRes; j++) {
            map[j] = new double[xRes];
            for (int i = 0; i < xRes; i++) {
                map[j][i] = 0;
            }
        }

        for (int octave = 1; octave <= numOctaves; octave++) {

            for (int j = 0; j < yRes; j++) {
                for (int i = 0; i < xRes; i++) {
                    double x = static_cast<double>(octave * i) / static_cast<double>(xRes) + octave;
                    double y = static_cast<double>(octave * j) / static_cast<double>(yRes) + octave;
                    map[j][i] += noise.getVal(x, y) / (1 << (octave - 1));
                }
            }

        }
        
        mapToRange(0.0, 255.0);
    }

    void mapToRange(double lowerBound, double upperBound) {
        double minVal;
        double maxVal;
        
        for (int j = 0; j < yRes; j++) {
            for (int i = 0; i < xRes; i++) {
                
                if (i == 0 && j == 0) {
                    minVal = map[j][i];
                    maxVal = map[j][i];
                }
                else {
                    if (map[j][i] < minVal) {
                        minVal = map[j][i];
                    }
                    if (map[j][i] > maxVal) {
                        maxVal = map[j][i];
                    }
                }
                
            }
        }

        for (int j = 0; j < yRes; j++) {
            for (int i = 0; i < xRes; i++) {
                // Map to range [0, 1]
                map[j][i] = (map[j][i] - minVal) / (maxVal - minVal);
                // Map to range [lowerBound, upperBound]
                map[j][i] = map[j][i] * (upperBound - lowerBound) + lowerBound;
            }
        }
    }

    double* operator[](int j) {
        return map[j];
    }

    int getXRes() {
        return xRes;
    }

    int getYRes() {
        return yRes;
    }

private:
    double** map;
    int xRes;
    int yRes;
};