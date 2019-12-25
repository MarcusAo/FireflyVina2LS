#ifndef KMEANS_H_
#define KMEANS_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

class Result
{
public:
    Result(int, int, int);
    int value1;
    int value2;
    int value3;
};

class Point
{

private:
    int pointId, clusterId;
    int dimensions;
    vector<double> values;
    double fit;

public:
    Point(int, string);

    Point(int, double);

    Point(int, vector<double>);

    int getDimensions();

    int getCluster();

    int getID();
    
    double getFit();

    void setCluster(int);

    double getVal(int);
};

class Cluster
{

private:
    int clusterId;
    vector<double> centroid;
    vector<Point> points;

public:
    Cluster(int, Point);

    void addPoint(Point);

    bool removePoint(int);

    int getId();

    Point getPoint(int);

    int getSize();

    double getCentroidByPos(int);

    void setCentroidByPos(int, double);
};

class KMeans
{
private:
    int K, iters, dimensions, total_points;
    vector<Cluster> clusters;

    int getNearestClusterId(Point);

public:
    KMeans(int, int);
    int a, b, c;
    void run(double[], int);

    void run(vector<double>[], int);
};

#endif