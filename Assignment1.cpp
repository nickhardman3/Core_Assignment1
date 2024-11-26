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
    double x; double y; //creates a point in the x and y axis
};

using utils::rand::randDouble; //using function from rng.cpp to create a random double value

double distancecalc(const point& a, const point&b)
{
    return std::sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)); //pythagoras theorem to work out the distance between two points
}

double wrapdistance(const point& a, const point&b)
{
    double dx=std::min(std::abs(a.x-b.x), 1.0-std::abs(a.x-b.x)); //computes the serial and wraparound distances and takes the smaller
    double dy=std::min(std::abs(a.y-b.y), 1.0-std::abs(a.y-b.y));

    return std::sqrt(dx*dx+dy*dy); //pythagoras
}

void randompoints(int N) //creating random points, comparing distances and working out averages
{
    std::vector<point> data; //creates a vector to store all locations in 2D
    data.reserve(N);

    utils::rand::seedRand(324255342);

    for (int i=0; i<N; ++i) //whilst i<number of locations
    {
        double x = utils::rand::randDouble(0.0, 1.0); //creates a random point between 0 and 1 for x and y axis
        double y = utils::rand::randDouble(0.0, 1.0);

        data.push_back({x, y}); //adds to the vector

    }

    std::vector<double> closestdist(N, std::numeric_limits<double>::max()); //stores closest and furthest distances
    std::vector<double> furthestdist(N, std::numeric_limits<double>::lowest());
    double closest_total=0.0, furthest_total=0.0; //variables which will add all distances to average

    double startTime = omp_get_wtime(); //starts timer

    for (int i=0; i<N; ++i)
    {
        for (int j=0; j<N; ++j) //eg, for every position i, gets compared to every position j
        {
            if (i==j) continue; // ignores when position i == j
            double dist = distancecalc(data[i], data[j]); //calc distance between points
            closestdist[i]=std::min(closestdist[i], dist); //takes the closest and adds to the closestdist vector
            furthestdist[i]=std::max(furthestdist[i], dist);//takes the furthest **
        }

        closest_total+=closestdist[i]; //adds the closest and furthest to the total variables
        furthest_total+=furthestdist[i];
    }

    double elapsedTime = omp_get_wtime() - startTime; //ends timer and works out the time elapsed

    std::string close_random = std::to_string(N) + "closest_distances_unparallelised_unoptimised.csv"; //creates csv starting with the number of locations used
    std::ofstream close(close_random);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl; //adds all closest distances to the csv
    }
    close.close();

    std::string far_random = std::to_string(N) + "furthest_distances_unparallelised_unoptimised.csv"; //** 
    std::ofstream far(far_random);
    for (int i=0; i<N; ++i)
    {
        far << furthestdist[i] << std::endl;
    }
    far.close();

    double close_avg = closest_total/N; //divides total by number of locations to work out averages
    double far_avg = furthest_total/N;

    std::cout << "Elapsed time: " << elapsedTime << "seconds" << std::endl;
    std::cout << "Average distance to the closest object: " << close_avg << std::endl;
    std::cout << "Average distance to the furthest object: " << far_avg << std::endl;
}

void readcsv(std::istream& in, std::vector<std::vector<double>>& data)
{
    std::string line; 
    while(std::getline(in, line)) //read the file by line
    {
        std::vector<double> row; //temporary vector 
        std::istringstream lineIn(line); //create string stream
        std::string cell;
        while(std::getline(lineIn, cell, ',')) //split the line by comma
        {
            row.push_back(std::stod(cell));  //convert string to double
        }
        data.push_back(row); //adds row to data vector
    }
}

void csvpoints(int N) //loading and processing csv files
{
    std::vector<point> data;
    std::string csv = std::to_string(N) + " locations.csv"; //selects the correct csv file given input of 100000 or 200000
    std::ifstream infile(csv);
    std::vector<std::vector<double>> raw_data;
    readcsv(infile, raw_data); //reads the csv file chosen
    infile.close();

    for (const auto& row: raw_data)
    {
        if (row.size () >=2) //ensures that the row has at least 2 elements
        {
            data.push_back({static_cast<double>(row[0]), static_cast<double>(row[1])}); //converts the first two elements of the row to doubles and create a point object
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

void wraprandompoints(int N) //wrap
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
            double dist = wrapdistance(data[i], data[j]); //uses wrapdistance instead of distancecalc
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

void wrapcsvpoints(int N) //no new comments
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

void parallelrandompoints(int N) //parallelisation
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

    #pragma omp parallel for reduction(+:closest_total,furthest_total) schedule(dynamic, 2) //parallelize the loop, reduction to combine totals across threads, dynamic scheduling for better load balancing
    for (int i=0; i<N; ++i)
    {
        double local_closest = std::numeric_limits<double>::max();
        double local_furthest = std::numeric_limits<double>::lowest();

        for (int j=0; j<N; ++j)
        {
            if (i==j) continue;
            double dist = distancecalc(data[i], data[j]);
            
            local_closest = std::min(local_closest, dist); //update local distances
            local_furthest = std::max(local_furthest, dist);
        }
        closestdist[i] = local_closest; //store the distances for i
        furthestdist[i] = local_furthest;

        closest_total+=closestdist[i]; 
        furthest_total+=furthestdist[i];
    }

    double elapsedTime = omp_get_wtime() - startTime;

    std::string close_random_para = std::to_string(N) + "closest_distances_parallelised_unoptimised.csv";
    std::ofstream close(close_random_para);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string far_random_para = std::to_string(N) + "furthest_distances_parallelised_unoptimised.csv";
    std::ofstream far(far_random_para);
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

void parallelcsvpoints(int N) //no new comments
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

    #pragma omp parallel for reduction(+:closest_total,furthest_total) schedule(dynamic, 2)
    for (int i=0; i<N; ++i)
    {
        double local_closest = std::numeric_limits<double>::max();
        double local_furthest = std::numeric_limits<double>::lowest();

        for (int j=0; j<N; ++j)
        {
            if (i==j) continue;
            double dist = distancecalc(data[i], data[j]);
            
            local_closest = std::min(local_closest, dist);
            local_furthest = std::max(local_furthest, dist);
        }
        closestdist[i] = local_closest;
        furthestdist[i] = local_furthest;

        closest_total+=closestdist[i];
        furthest_total+=furthestdist[i];
    }

    double elapsedTime = omp_get_wtime() - startTime;

    std::string csv_close_para = std::to_string(N) + "csv_closest_distances_parallelised_unoptimised.csv";
    std::ofstream close(csv_close_para);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string csv_far_para = std::to_string(N) + "csv_furthest_distances_parallelised_unoptimised.csv";
    std::ofstream far(csv_far_para);
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

void parallelwraprandompoints(int N) //no new comments
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

    #pragma omp parallel for reduction(+:closest_total,furthest_total) schedule(dynamic, 2)
    for (int i=0; i<N; ++i)
    {
        double local_closest = std::numeric_limits<double>::max();
        double local_furthest = std::numeric_limits<double>::lowest();

        for (int j=0; j<N; ++j)
        {
            if (i==j) continue;
            double dist = wrapdistance(data[i], data[j]);
            
            local_closest = std::min(local_closest, dist);
            local_furthest = std::max(local_furthest, dist);
        }
        closestdist[i] = local_closest;
        furthestdist[i] = local_furthest;

        closest_total+=closestdist[i];
        furthest_total+=furthestdist[i];
    }

    double elapsedTime = omp_get_wtime() - startTime;

    std::string close_random_para = std::to_string(N) + "closest_distances_parallelised_wrap.csv";
    std::ofstream close(close_random_para);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string far_random_para = std::to_string(N) + "furthest_distances_parallelised_wrap.csv";
    std::ofstream far(far_random_para);
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

void parallelwrapcsvpoints(int N) //no new comments
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

    #pragma omp parallel for reduction(+:closest_total,furthest_total) schedule(dynamic, 2)
    for (int i=0; i<N; ++i)
    {
        double local_closest = std::numeric_limits<double>::max();
        double local_furthest = std::numeric_limits<double>::lowest();

        for (int j=0; j<N; ++j)
        {
            if (i==j) continue;
            double dist = wrapdistance(data[i], data[j]);
            
            local_closest = std::min(local_closest, dist);
            local_furthest = std::max(local_furthest, dist);
        }
        closestdist[i] = local_closest;
        furthestdist[i] = local_furthest;

        closest_total+=closestdist[i];
        furthest_total+=furthestdist[i];
    }

    double elapsedTime = omp_get_wtime() - startTime;

    std::string csv_close = std::to_string(N) + "csv_closest_distances_parallelised_wrap.csv";
    std::ofstream close(csv_close);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string csv_far = std::to_string(N) + "csv_furthest_distances_parallelised_wrap.csv";
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

void optimisedrandompoints(int N) //optimising 
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
        for (int j=i+1; j<N; ++j) //j=i+1 means that distance calculated once as eg, distance wont be calculated from 2 to 1
        {
            double dist = distancecalc(data[i], data[j]);

            closestdist[i]=std::min(closestdist[i], dist); //updates distances for both i and j based on calculated distances
            furthestdist[i]=std::max(furthestdist[i], dist);

            closestdist[j]=std::min(closestdist[j], dist);
            furthestdist[j]=std::max(furthestdist[j], dist);            
        }

        closest_total+=closestdist[i];
        furthest_total+=furthestdist[i];
    }

    double elapsedTime = omp_get_wtime() - startTime;

    std::string close_random = std::to_string(N) + "closest_distances_unparallelised_optimised.csv";
    std::ofstream close(close_random);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string far_random = std::to_string(N) + "furthest_distances_unparallelised_optimised.csv";
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

void optimisedcsvpoints(int N) //no new comments
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
        for (int j=i+1; j<N; ++j)
        {
            double dist = distancecalc(data[i], data[j]);

            closestdist[i]=std::min(closestdist[i], dist);
            furthestdist[i]=std::max(furthestdist[i], dist);

            closestdist[j]=std::min(closestdist[j], dist);
            furthestdist[j]=std::max(furthestdist[j], dist);            
        }

        closest_total+=closestdist[i];
        furthest_total+=furthestdist[i];
    }

    double elapsedTime = omp_get_wtime() - startTime;

    std::string csv_close = std::to_string(N) + "csv_closest_distances_unparallelised_optimised.csv";
    std::ofstream close(csv_close);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string csv_far = std::to_string(N) + "csv_furthest_distances_unparallelised_optimised.csv";
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

void optimisedwraprandompoints(int N) //no new comments
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
        for (int j=i+1; j<N; ++j)
        {
            double dist = distancecalc(data[i], data[j]);

            closestdist[i]=std::min(closestdist[i], dist);
            furthestdist[i]=std::max(furthestdist[i], dist);

            closestdist[j]=std::min(closestdist[j], dist);
            furthestdist[j]=std::max(furthestdist[j], dist);            
        }

        closest_total+=closestdist[i];
        furthest_total+=furthestdist[i];
    }

    double elapsedTime = omp_get_wtime() - startTime;

    std::string close_random = std::to_string(N) + "closest_distances_unparallelised_wrap_optimised.csv";
    std::ofstream close(close_random);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string far_random = std::to_string(N) + "furthest_distances_unparallelised_wrap_optimised.csv";
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

void optimisedwrapcsvpoints(int N) //no new comments
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
        for (int j=i+1; j<N; ++j)
        {
            double dist = distancecalc(data[i], data[j]);

            closestdist[i]=std::min(closestdist[i], dist);
            furthestdist[i]=std::max(furthestdist[i], dist);

            closestdist[j]=std::min(closestdist[j], dist);
            furthestdist[j]=std::max(furthestdist[j], dist);            
        }

        closest_total+=closestdist[i];
        furthest_total+=furthestdist[i];
    }

    double elapsedTime = omp_get_wtime() - startTime;

    std::string csv_close = std::to_string(N) + "csv_closest_distances_unparallelised_wrap_optimised.csv";
    std::ofstream close(csv_close);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string csv_far = std::to_string(N) + "csv_furthest_distances_unparallelised_wrap_optimised.csv";
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

void optimisedparallelrandompoints(int N) //optimising parallelisation
{
    std::vector<point> data;
    data.reserve(N);

    utils::rand::seedRand(324255342);

    for (int i = 0; i < N; ++i) 
    {
        double x = utils::rand::randDouble(0.0, 1.0);
        double y = utils::rand::randDouble(0.0, 1.0);
        data.push_back({x, y});
    }

    std::vector<double> closestdist(N, std::numeric_limits<double>::max());
    std::vector<double> furthestdist(N, std::numeric_limits<double>::lowest());

    double closest_total = 0.0, furthest_total = 0.0;

    double startTime = omp_get_wtime();

    #pragma omp parallel
    {
        std::vector<double> local_closest(N, std::numeric_limits<double>::max()); //local storage for intermediate close & far distances
        std::vector<double> local_furthest(N, std::numeric_limits<double>::lowest());

        #pragma omp for reduction(+:closest_total, furthest_total) schedule(dynamic, 2)
        for (int i = 0; i < N; ++i) 
        {
            for (int j = i + 1; j < N; ++j) 
            {
                double dist = distancecalc(data[i], data[j]);

                local_closest[i] = std::min(local_closest[i], dist);
                local_furthest[i] = std::max(local_furthest[i], dist);

                local_closest[j] = std::min(local_closest[j], dist);
                local_furthest[j] = std::max(local_furthest[j], dist);
            }
        }

        #pragma omp critical //critical region used to update the global close&far distances
        {
            for (int i = 0; i < N; ++i) 
            {
                closestdist[i] = std::min(closestdist[i], local_closest[i]);
                furthestdist[i] = std::max(furthestdist[i], local_furthest[i]);
            }
        }
    }

    for (int i = 0; i < N; ++i) 
    {
        closest_total += closestdist[i];
        furthest_total += furthestdist[i];
    }

    double elapsedTime = omp_get_wtime() - startTime;

    std::string close_random_para = std::to_string(N) + "closest_distances_parallelised_optimised.csv";
    std::ofstream close(close_random_para);
    for (int i = 0; i < N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string far_random_para = std::to_string(N) + "furthest_distances_parallelised_optimised.csv";
    std::ofstream far(far_random_para);
    for (int i = 0; i < N; ++i)
    {
        far << furthestdist[i] << std::endl;
    }
    far.close();

    double close_avg = closest_total / N;
    double far_avg = furthest_total / N;

    std::cout << "Elapsed time: " << elapsedTime << " seconds" << std::endl;
    std::cout << "Average distance to the closest object: " << close_avg << std::endl;
    std::cout << "Average distance to the furthest object: " << far_avg << std::endl;

}

void optimisedparallelcsvpoints(int N) //no new comments
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

    #pragma omp parallel
    {
        std::vector<double> local_closest(N, std::numeric_limits<double>::max());
        std::vector<double> local_furthest(N, std::numeric_limits<double>::lowest());

        #pragma omp for reduction(+:closest_total, furthest_total) schedule(dynamic, 2)
        for (int i = 0; i < N; ++i) 
        {
            for (int j = i + 1; j < N; ++j) 
            {
                double dist = distancecalc(data[i], data[j]);

                local_closest[i] = std::min(local_closest[i], dist);
                local_furthest[i] = std::max(local_furthest[i], dist);

                local_closest[j] = std::min(local_closest[j], dist);
                local_furthest[j] = std::max(local_furthest[j], dist);
            }
        }

        #pragma omp critical
        {
            for (int i = 0; i < N; ++i) 
            {
                closestdist[i] = std::min(closestdist[i], local_closest[i]);
                furthestdist[i] = std::max(furthestdist[i], local_furthest[i]);
            }
        }
    }

    for (int i = 0; i < N; ++i) 
    {
        closest_total += closestdist[i];
        furthest_total += furthestdist[i];
    }

    double elapsedTime = omp_get_wtime() - startTime;

    std::string csv_close_para = std::to_string(N) + "csv_closest_distances_parallelised_unoptimised.csv";
    std::ofstream close(csv_close_para);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string csv_far_para = std::to_string(N) + "csv_furthest_distances_parallelised_unoptimised.csv";
    std::ofstream far(csv_far_para);
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

void optimisedparallelwraprandompoints(int N) //no new comments
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

    #pragma omp parallel
    {
        std::vector<double> local_closest(N, std::numeric_limits<double>::max());
        std::vector<double> local_furthest(N, std::numeric_limits<double>::lowest());

        #pragma omp for reduction(+:closest_total, furthest_total) schedule(dynamic, 2)
        for (int i = 0; i < N; ++i) 
        {
            for (int j = i + 1; j < N; ++j) 
            {
                double dist = wrapdistance(data[i], data[j]);

                local_closest[i] = std::min(local_closest[i], dist);
                local_furthest[i] = std::max(local_furthest[i], dist);

                local_closest[j] = std::min(local_closest[j], dist);
                local_furthest[j] = std::max(local_furthest[j], dist);
            }
        }

        #pragma omp critical
        {
            for (int i = 0; i < N; ++i) 
            {
                closestdist[i] = std::min(closestdist[i], local_closest[i]);
                furthestdist[i] = std::max(furthestdist[i], local_furthest[i]);
            }
        }
    }

    for (int i = 0; i < N; ++i) 
    {
        closest_total += closestdist[i];
        furthest_total += furthestdist[i];
    }


    double elapsedTime = omp_get_wtime() - startTime;

    std::string close_random_para = std::to_string(N) + "closest_distances_parallelised_wrap.csv";
    std::ofstream close(close_random_para);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string far_random_para = std::to_string(N) + "furthest_distances_parallelised_wrap.csv";
    std::ofstream far(far_random_para);
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

void optimisedparallelwrapcsvpoints(int N) //no new comments
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

    #pragma omp parallel
    {
        std::vector<double> local_closest(N, std::numeric_limits<double>::max());
        std::vector<double> local_furthest(N, std::numeric_limits<double>::lowest());

        #pragma omp for reduction(+:closest_total, furthest_total) schedule(dynamic, 2)
        for (int i = 0; i < N; ++i) 
        {
            for (int j = i + 1; j < N; ++j) 
            {
                double dist = wrapdistance(data[i], data[j]);

                local_closest[i] = std::min(local_closest[i], dist);
                local_furthest[i] = std::max(local_furthest[i], dist);

                local_closest[j] = std::min(local_closest[j], dist);
                local_furthest[j] = std::max(local_furthest[j], dist);
            }
        }

        #pragma omp critical
        {
            for (int i = 0; i < N; ++i) 
            {
                closestdist[i] = std::min(closestdist[i], local_closest[i]);
                furthestdist[i] = std::max(furthestdist[i], local_furthest[i]);
            }
        }
    }

    for (int i = 0; i < N; ++i) 
    {
        closest_total += closestdist[i];
        furthest_total += furthestdist[i];
    }

    double elapsedTime = omp_get_wtime() - startTime;

    std::string csv_close = std::to_string(N) + "csv_closest_distances_parallelised_wrap.csv";
    std::ofstream close(csv_close);
    for (int i=0; i<N; ++i)
    {
        close << closestdist[i] << std::endl;
    }
    close.close();

    std::string csv_far = std::to_string(N) + "csv_furthest_distances_parallelised_wrap.csv";
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

    if(choice==1) //choices leading to a function of choice being ran
    {
        std::cout << "100000 or 200000 object locations?";
        int N;
        std::cin >> N;
        if (N==100000 || N==200000) //choose how many locations
        {
            std::cout << "Choose calculation (1 = Serial, 2 = Wraparound)";
            int calc;
            std::cin >> calc;

            if (calc == 1 || calc == 2) //serial or wrap
            {
                std::cout << "Optimised? (1 = No, 2 = Yes)";
                int opt;
                std::cin >> opt;

                if (opt == 1 || opt == 2) //optimised or not
                {
                    std::cout << "Parallelised? (1 = No, 2 = Yes)";
                    int para;
                    std::cin >> para;

                    if (para == 1) //parallelised or not
                    {
                        if (opt == 1)
                        {
                            if (calc == 1)
                            {
                                randompoints(N);
                            }
                            else if (calc == 2)
                            {
                                wraprandompoints(N);
                            }
                        }
                        else if (opt == 2)
                        {
                            if (calc == 1)
                            {
                                optimisedrandompoints(N);
                            }
                            else if (calc == 2)
                            {
                                optimisedwraprandompoints(N);
                            }
                        }
                    }
                    else if (para == 2)
                    {
                        std::cout << "How many threads? (1-8)"; //gives a choice of how many threads to use in parallelisation
                        int threads;
                        std::cin >> threads;
                        omp_set_num_threads(threads);

                        if (opt == 1)
                        {
                            if (calc == 1)
                            {
                                parallelrandompoints(N);
                            }
                            else if (calc == 2)
                            {
                                parallelwraprandompoints(N);
                            }
                        }
                        else if (opt == 2)
                        {
                            if (calc == 1)
                            {
                                optimisedparallelrandompoints(N);
                            }
                            else if (calc == 2)
                            {
                                optimisedparallelwraprandompoints(N);
                            }
                        }
                    } 
                    else
                    {
                        std::cerr << "Invalid choice";
                    }                               
                }
                else
                {
                    std::cerr << "Invalid choice";
                }
            }
            else
            {
                std::cerr << "Invalid choice";
            }
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
            std::cout << "Choose calculation (1 = Serial, 2 = Wraparound)";
            int calc;
            std::cin >> calc;

            if (calc == 1 || calc == 2)
            {
                std::cout << "Optimised? (1 = No, 2 = Yes)";
                int opt;
                std::cin >> opt;

                if (opt == 1 || opt == 2)
                {
                    std::cout << "Parallelised? (1 = No, 2 = Yes)";
                    int para;
                    std::cin >> para;

                    if (para == 1)
                    {
                        if (opt == 1)
                        {
                            if (calc == 1)
                            {
                                csvpoints(N);
                            }
                            else if (calc == 2)
                            {
                                wrapcsvpoints(N);
                            }
                        }
                        else if (opt == 2)
                        {
                            if (calc == 1)
                            {
                                optimisedcsvpoints(N);
                            }
                            else if (calc == 2)
                            {
                                optimisedwrapcsvpoints(N);
                            }
                        }
                    }
                    else if (para == 2)
                    {
                        std::cout << "How many threads? (1-8)";
                        int threads;
                        std::cin >> threads;
                        omp_set_num_threads(threads);

                        if (opt == 1)
                        {
                            if (calc == 1)
                            {
                                parallelcsvpoints(N);
                            }
                            else if (calc == 2)
                            {
                                parallelwrapcsvpoints(N);
                            }
                        }
                        else if (opt == 2)
                        {
                            if (calc == 1)
                            {
                                optimisedparallelcsvpoints(N);
                            }
                            else if (calc == 2)
                            {
                                optimisedparallelwrapcsvpoints(N);
                            }
                        }
                    } 
                    else
                    {
                        std::cerr << "Invalid choice";
                    }                               
                }
                else
                {
                    std::cerr << "Invalid choice";
                }
            }
            else
            {
                std::cerr << "Invalid choice";
            }
        }
        else
        {
            std::cerr << "Invalid choice";
        }        
    }
    else
    {
        std::cerr << "Invalid choice";
    }
    
    return EXIT_SUCCESS;
}
