#pragma once
#include <cstdlib>
#include <cmath>

class PerlinNoise {
public:
    PerlinNoise(int p) {
        period = p;
        grid = new int* [period];
        for (int i = 0; i < period; i++) {
            grid[i] = new int[period];
            for (int j = 0; j < period; j++) {
                grid[i][j] = rand() % 8;
            }
        }
    }

    double getVal(double x, double y) {
        x = fmod(x, static_cast<double>(period));
        y = fmod(y, static_cast<double>(period));
        int minX = static_cast<int>(floor(x));
        int minY = static_cast<int>(floor(y));
        int maxX = minX + 1;
        int maxY = minY + 1;

        double fracX = x - minX;
        double fracY = y - minY;

        double d00[2] = { fracX, fracY };
        double d10[2] = { fracX - 1.0, fracY };
        double d01[2] = { fracX, fracY - 1.0 };
        double d11[2] = { fracX - 1.0, fracY - 1.0 };

        double n00 = dot(gradients[grid[minX][minY]], d00);
        double n10 = dot(gradients[grid[maxX][minY]], d10);
        double n01 = dot(gradients[grid[minX][maxY]], d01);
        double n11 = dot(gradients[grid[maxX][maxY]], d11);

        double fadeX = fade(fracX);
        double fadeY = fade(fracY);

        double nx0 = n00 * (1 - fadeX) + n10 * fadeX;
        double nx1 = n01 * (1 - fadeX) + n11 * fadeX;

        double nxy = nx0 * (1 - fadeY) + nx1 * fadeY;

        return nxy;
    }

private:
    double dot(double vec0[2], double vec1[2]) {
        return vec0[0] * vec1[0] + vec0[1] * vec1[1];
    }

    double fade(double t) {
        return 6 * pow(t, 5.0) - 15 * pow(t, 4.0) + 10 * pow(t, 3.0);
    }

    //double gradients[4][2] = { {1.0, 1.0}, {1.0, -1.0}, {-1.0, -1.0}, {-1.0, 1.0} };
    double gradients[8][2] = { {sqrt(2.0), 0}, {0, sqrt(2.0)}, {-sqrt(2.0), 0}, {0, -sqrt(2.0)}, {1.0, 1.0}, {1.0, -1.0}, {-1.0, -1.0}, {-1.0, 1.0} };
    int** grid;
    int period;
};