#pragma once
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "mpi.h"

#include "ParallelHeightMap.h"

// Different erosion pattern than serial?

class ParallelDrop {
public:
    ParallelDrop(ParallelHeightMap* m, double positionX, double positionY) {
        map = m;
        position[0] = positionX;
        position[1] = positionY;
        direction[0] = 0.0;
        direction[1] = 0.0;
        speed = 0.0;
        water = 1.0;
        sediment = 0.0;
        capacity = 0.0;
        numIters = 0;
    }

    bool iterate() {
        if (position[0] >= map->getXRes() - 1.0 || position[1] >= map->getYRes() - 1.0 ||
            position[0] < 0 || position[1] < 0 || numIters >= MAXIMUM_PATH) {
            return false;
        }

        numIters++;

        MPI_Win win = map->getWin();

        int x0 = static_cast<int>(position[0]);
        int y0 = static_cast<int>(position[1]);
        int x1 = x0 + 1;
        int y1 = y0 + 1;

        // Get surrounding heights from shared memory
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
        double h00 = (*map)[y0][x0];
        double h10 = (*map)[y0][x1];
        double h01 = (*map)[y1][x0];
        double h11 = (*map)[y1][x1];
        MPI_Win_unlock(0, win);
        
        // Get fractional distance from x0 and y0
        double fracX = position[0] - x0;
        double fracY = position[1] - y0;

        // Height changes of nearest points
        double deltaX0 = h10 - h00;
        double deltaX1 = h11 - h01;
        double delta0Y = h01 - h00;
        double delta1Y = h11 - h10;

        // Negative gradient at current position
        double gcur[2] = { -mix(deltaX0, deltaX1, fracY), -mix(delta0Y, delta1Y, fracX) };

        // Get new direction
        direction[0] = mix(gcur[0], direction[0], INERTIA);
        direction[1] = mix(gcur[1], direction[1], INERTIA);

        // Choose random direction if magnitude of new direction is 0.0
        if (direction[0] == 0.0 && direction[1] == 0.0) {
            direction[0] = randomDouble();
            direction[1] = randomDouble();
        }

        normalize(direction);

        // Interpolate height at old position
        double oldHeight = mix(mix(h00, h10, fracX), mix(h01, h11, fracX), fracY);

        // Get new position
        double newPosition[2] = { position[0] + direction[0], position[1] + direction[1] };

        // Get new surrounding points
        int nx0 = static_cast<int>(newPosition[0]);
        int ny0 = static_cast<int>(newPosition[1]);
        int nx1 = x0 + 1;
        int ny1 = y0 + 1;

        // Get new surrounding heights from shared memory
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
        double nh00 = (*map)[ny0][nx0];
        double nh10 = (*map)[ny0][nx1];
        double nh01 = (*map)[ny1][nx0];
        double nh11 = (*map)[ny1][nx1];
        MPI_Win_unlock(0, win);
        
        // Interpolate height at new position
        double newHeight = mix(mix(nh00, nh10, fracX), mix(nh01, nh11, fracX), fracY);

        // Get change in height
        double heightDiff = newHeight - oldHeight;
        
        // Track amounts to deposit
        double dep00 = 0.0;
        double dep10 = 0.0;
        double dep01 = 0.0;
        double dep11 = 0.0;

        if (heightDiff > 0.0) {
            // Deposit heightDiff sediment if possible, otherwise deposit all sediment
            double amount = std::min(heightDiff, sediment);

            dep00 += (1 - fracX) * (1 - fracY) * amount;
            dep10 += fracX * (1 - fracY) * amount;
            dep01 += (1 - fracX) * fracY * amount;
            dep11 += fracX * fracY * amount;
            
            sediment -= amount;

            // Reset speed to zero
            speed = 0.0;

        }
        else {
            // Adjust capacity if going downhill
            capacity = std::max(-heightDiff, MINIMUM_SLOPE) * speed * water * CAPACITY;

            // Get new speed
            speed = sqrt(speed * speed + heightDiff * GRAVITY);
        }

        if (sediment > capacity) {
            // Drop (sediment - capacity) * DEPOSITION_SPEED at old position (this value is always less than sediment)
            double amount = (sediment - capacity) * DEPOSITION_SPEED;

            dep00 += (1 - fracX) * (1 - fracY) * amount;
            dep10 += fracX * (1 - fracY) * amount;
            dep01 += (1 - fracX) * fracY * amount;
            dep11 += fracX * fracY * amount;

            sediment -= amount;
        }
        else if (sediment < capacity) {
            // Take std::min((capacity - sediment) * EROSION_SPEED, -heightDiff) from points within EROSION_RADIUS of old position
            double amount = (capacity - sediment) * EROSION_SPEED;
            if (heightDiff < 0) {
                amount = std::min(amount, -heightDiff);
            }

            int minX = std::max(static_cast<int>(position[0] - EROSION_RADIUS), 0);
            int maxX = std::min(static_cast<int>(position[0] + EROSION_RADIUS), map->getXRes() - 1);
            int minY = std::max(static_cast<int>(position[1] - EROSION_RADIUS), 0);
            int maxY = std::min(static_cast<int>(position[1] + EROSION_RADIUS), map->getYRes() - 1);

            double total = 0.0;
            double** weights = new double* [maxY - minY + 1];
            for (int j = 0; j < maxY - minY + 1; j++) {
                weights[j] = new double[maxX - minX + 1];
                for (int i = 0; i < maxX - minX + 1; i++) {
                    double xDiff = static_cast<double>(i + minX) - position[0];
                    double yDiff = static_cast<double>(j + minY) - position[1];
                    weights[j][i] = std::max(0.0, EROSION_RADIUS - sqrt(xDiff * xDiff + yDiff * yDiff));
                    total += weights[j][i];
                }
            }

            // Put changes into shared memory
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
            for (int j = 0; j < maxY - minY + 1; j++) {
                for (int i = 0; i < maxX - minX + 1; i++) {
                    (*map)[j + minY][i + minX] -= amount * weights[j][i] / total;
                }
            }
            
            (*map)[y0][x0] += dep00;
            (*map)[y0][x1] += dep10;
            (*map)[y1][x0] += dep01;
            (*map)[y1][x1] += dep11;
            MPI_Win_unlock(0, win);

            sediment += amount;
        }
        
        // Get new water amount
        water *= 1.0 - EVAPORATION_SPEED;

        // Set new position
        position[0] = newPosition[0];
        position[1] = newPosition[1];

        return true;
    }

    double getPosX() {
        return position[0];
    }

    double getPosY() {
        return position[1];
    }

private:
    double mix(double v0, double v1, double frac) {
        return v0 * (1 - frac) + v1 * frac;
    }

    void normalize(double v[2]) {
        double magnitude = sqrt(v[0] * v[0] + v[1] * v[1]);
        v[0] /= magnitude;
        v[1] /= magnitude;
    }

    // Random double between -1.0 and 1.0
    double randomDouble() {
        return static_cast<double>(rand() / RAND_MAX) * 2.0 - 1.0;
    }

    static const double INERTIA;
    static const double CAPACITY;
    static const double DEPOSITION_SPEED;
    static const double EROSION_SPEED;
    static const double EVAPORATION_SPEED;
    static const double EROSION_RADIUS;
    static const double MINIMUM_SLOPE;
    static const int    MAXIMUM_PATH;
    static const double GRAVITY;

    ParallelHeightMap* map;
    double position[2];
    double direction[2];
    double speed;
    double water;
    double sediment;
    double capacity;
    int numIters;
};


const double ParallelDrop::INERTIA = 0.1;
const double ParallelDrop::CAPACITY = 8.0;
const double ParallelDrop::DEPOSITION_SPEED = 0.1;
const double ParallelDrop::EROSION_SPEED = 0.1;
const double ParallelDrop::EVAPORATION_SPEED = 0.05;
const double ParallelDrop::EROSION_RADIUS = 5.0;
const double ParallelDrop::MINIMUM_SLOPE = 0.01;
const int    ParallelDrop::MAXIMUM_PATH = 30;
const double ParallelDrop::GRAVITY = -10.0;