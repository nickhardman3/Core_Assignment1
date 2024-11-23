#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include "rng.h"
#include <cmath>
#include <limits>

struct point 
{
    double x; double y;
};

using utils::rand::randDouble;

double distancecalc(const point& a, const point&b)
{
    return std::sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}

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

//    std::ofstream out("particle_points_unoptimised_unparralisd.csv");

//    for (int i=0; i<data.size(); ++i)
//    {
//        out << data[i].x << "," << data[i].y << std::endl;
//    }
//    out.close();

    std::vector<double> closestdist(N, std::numeric_limits<double>::max());
    std::vector<double> furthestdist(N, std::numeric_limits<double>::min());
    double closest_total=0.0, furthest_total=0.0;

    for (int i=0; i<N; ++i)
    {
        for (int j=0; j<N; ++j)
        {
            if (i==j) continue;
            double dist = distancecalc(data[i], data[j]);
            closestdist[i]=std::min(closestdist[i], dist);
            furthestdist[i]=std::max(furthestdist[i], dist);
        }

        closest_total+=closestdist[i];
        furthest_total+=furthestdist[i];
    }

    std::ofstream close("closest_distances_unparallelised_unoptimised.csv");
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::ofstream far("furthest_distances_unparallelised_unoptimised.csv");
    for (int i=0; i<N; ++i)
    {
        far << furthestdist[i] << std::endl;
    }
    far.close();

    double close_avg = closest_total/N;
    double far_avg = furthest_total/N;

    std::cout << "Average distance to the closest object: " << close_avg << std::endl;
    std::cout << "Average distance to the furthest object: " << far_avg << std::endl;


    return EXIT_SUCCESS;
}