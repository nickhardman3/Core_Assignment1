#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include "rng.h"
#include <cmath>
#include <limits>
#include <omp.h>
#include <sstream>

struct point 
{
    double x; double y;
};

using utils::rand::randDouble;

double distancecalc(const point& a, const point&b)
{
    return std::sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}

double wrapdistance(const point& a, const point&b)
{
    double dx=std::min(std::abs(a.x-b.x), 1.0-std::abs(a.x-b.x));
    double dy=std::min(std::abs(a.y-b.y), 1.0-std::abs(a.y-b.y));

    return std::sqrt(dx*dx+dy*dy);
}

void randompoints(int N) 
{
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
    std::vector<double> furthestdist(N, std::numeric_limits<double>::lowest());
    double closest_total=0.0, furthest_total=0.0;

    double startTime = omp_get_wtime();

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

    double elapsedTime = omp_get_wtime() - startTime;

    std::string close_random = std::to_string(N) + "closest_distances_unparallelised_unoptimised.csv";
    std::ofstream close(close_random);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string far_random = std::to_string(N) + "furthest_distances_unparallelised_unoptimised.csv";
    std::ofstream far(far_random);
    for (int i=0; i<N; ++i)
    {
        far << furthestdist[i] << std::endl;
    }
    far.close();

    double close_avg = closest_total/N;
    double far_avg = furthest_total/N;

    std::cout << "Elapsed time: " << elapsedTime << "seconds" << std::endl;
    std::cout << "Average distance to the closest object: " << close_avg << std::endl;
    std::cout << "Average distance to the furthest object: " << far_avg << std::endl;
}

void readcsv(std::istream& in, std::vector<std::vector<double>>& data)
{
    std::string line;
    while(std::getline(in, line))
    {
        std::vector<double> row;
        std::istringstream lineIn(line);
        std::string cell;
        while(std::getline(lineIn, cell, ','))
        {
            row.push_back(std::stod(cell));
        }
        data.push_back(row);
    }
}

void csvpoints(int N)
{
    std::vector<point> data;
    std::string csv = std::to_string(N) + " locations.csv";
    std::ifstream infile(csv);
    std::vector<std::vector<double>> raw_data;
    readcsv(infile, raw_data);
    infile.close();

    for (const auto& row: raw_data)
    {
        if (row.size () >=2)
        {
            data.push_back({static_cast<double>(row[0]), static_cast<double>(row[1])});
        }
    }

    std::vector<double> closestdist(N, std::numeric_limits<double>::max());
    std::vector<double> furthestdist(N, std::numeric_limits<double>::lowest());
    double closest_total=0.0, furthest_total=0.0;

    double startTime = omp_get_wtime();

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

    double elapsedTime = omp_get_wtime() - startTime;

    std::string csv_close = std::to_string(N) + "csv_closest_distances_unparallelised_unoptimised.csv";
    std::ofstream close(csv_close);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string csv_far = std::to_string(N) + "csv_furthest_distances_unparallelised_unoptimised.csv";
    std::ofstream far(csv_far);
    for (int i=0; i<N; ++i)
    {
        far << furthestdist[i] << std::endl;
    }
    far.close(); 

    double close_avg = closest_total/N;
    double far_avg = furthest_total/N;

    std::cout << "Elapsed time: " << elapsedTime << "seconds" << std::endl;
    std::cout << "Average distance to the closest object: " << close_avg << std::endl;
    std::cout << "Average distance to the furthest object: " << far_avg << std::endl;
}

void wraprandompoints(int N)
{
    std::vector<point> data;
    data.reserve(N);

    utils::rand::seedRand(324255342);

    for (int i=0; i<N; ++i) 
    {
        double x = utils::rand::randDouble(0.0, 1.0);
        double y = utils::rand::randDouble(0.0, 1.0);

        data.push_back({x, y});

    }

    std::vector<double> closestdist(N, std::numeric_limits<double>::max());
    std::vector<double> furthestdist(N, std::numeric_limits<double>::lowest());
    double closest_total=0.0, furthest_total=0.0;

    double startTime = omp_get_wtime();

    for (int i=0; i<N; ++i)
    {
        for (int j=0; j<N; ++j)
        {
            if (i==j) continue;
            double dist = wrapdistance(data[i], data[j]);
            closestdist[i]=std::min(closestdist[i], dist);
            furthestdist[i]=std::max(furthestdist[i], dist);
        }

        closest_total+=closestdist[i];
        furthest_total+=furthestdist[i];
    }

    double elapsedTime = omp_get_wtime() - startTime;

    std::string close_random = std::to_string(N) + "closest_distances_unparallelised_wrap.csv";
    std::ofstream close(close_random);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string far_random = std::to_string(N) + "furthest_distances_unparallelised_wrap.csv";
    std::ofstream far(far_random);
    for (int i=0; i<N; ++i)
    {
        far << furthestdist[i] << std::endl;
    }
    far.close();

    double close_avg = closest_total/N;
    double far_avg = furthest_total/N;

    std::cout << "Elapsed time: " << elapsedTime << "seconds" << std::endl;
    std::cout << "Average distance to the closest object: " << close_avg << std::endl;
    std::cout << "Average distance to the furthest object: " << far_avg << std::endl;

}

void wrapcsvpoints(int N)
{
    std::vector<point> data;
    std::string csv = std::to_string(N) + " locations.csv";
    std::ifstream infile(csv);
    std::vector<std::vector<double>> raw_data;
    readcsv(infile, raw_data);
    infile.close();

    for (const auto& row: raw_data)
    {
        if (row.size () >=2)
        {
            data.push_back({static_cast<double>(row[0]), static_cast<double>(row[1])});
        }
    }

    std::vector<double> closestdist(N, std::numeric_limits<double>::max());
    std::vector<double> furthestdist(N, std::numeric_limits<double>::lowest());
    double closest_total=0.0, furthest_total=0.0;

    double startTime = omp_get_wtime();

    for (int i=0; i<N; ++i)
    {
        for (int j=0; j<N; ++j)
        {
            if (i==j) continue;
            double dist = wrapdistance(data[i], data[j]);
            closestdist[i]=std::min(closestdist[i], dist);
            furthestdist[i]=std::max(furthestdist[i], dist);
        }

        closest_total+=closestdist[i];
        furthest_total+=furthestdist[i];
    }

    double elapsedTime = omp_get_wtime() - startTime;

    std::string csv_close = std::to_string(N) + "csv_closest_distances_unparallelised_wrap.csv";
    std::ofstream close(csv_close);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string csv_far = std::to_string(N) + "csv_furthest_distances_unparallelised_wrap.csv";
    std::ofstream far(csv_far);
    for (int i=0; i<N; ++i)
    {
        far << furthestdist[i] << std::endl;
    }
    far.close(); 

    double close_avg = closest_total/N;
    double far_avg = furthest_total/N;

    std::cout << "Elapsed time: " << elapsedTime << "seconds" << std::endl;
    std::cout << "Average distance to the closest object: " << close_avg << std::endl;
    std::cout << "Average distance to the furthest object: " << far_avg << std::endl;
}

int main()
{
    std::cout << "Choose input type (1 = Random object locations, 2 = CSV object locations)";
    int choice;
    std::cin >> choice;

    if(choice==1)
    {
        std::cout << "100000 or 200000 object locations?";
        int N;
        std::cin >> N;
        if (N==100000 || N==200000)
        {
            randompoints(N);
        }
        else 
        {
            std::cerr << "Invalid number";
        }
    }
    else if (choice == 2)
    {
        std::cout << "100000 or 200000 object locations?";
        int N;
        std::cin >> N;
        if (N==100000 || N==200000)
        {
            csvpoints(N);
        }
        else 
        {
            std::cerr << "Invalid number";
        }
    }

    else
    {
        std::cerr << "Invalid choice";
    }
    
    return EXIT_SUCCESS;
}