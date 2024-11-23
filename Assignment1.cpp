#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include "rng.h"

struct point 
{
    double x; double y;
};

using utils::rand::randDouble;

int main() 
{
    const int N = 100;
    std::vector<point> data;
    data.reserve(N);

    utils::rand::seedRand(324255342);

    for (int i=0; i<N; ++i) 
    {
        double x = utils::rand::randDouble(0.0, 1.0);
        double y = utils::rand::randDouble(0.0, 1.0);

        data.push_back({x, y});

    }

    std::ofstream out("particle_points_unoptimised_unparralisd.csv");

    for (int i=0; i<data.size(); ++i)
    {
        out << data[i].x << "," << data[i].y << std::endl;
    }
    out.close();

    return EXIT_SUCCESS;
}