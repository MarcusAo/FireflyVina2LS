#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "kmeans.h"
using namespace std;

Result::Result(int a, int b, int c)
{
    value1 = a;
    value2 = b;
    value3 = c;
}

Point::Point(int id, string line)
{
    dimensions = 0;
    pointId = id;
    stringstream is(line);
    double val;
    while (is >> val)
    {
        values.push_back(val);
        dimensions++;
    }
    clusterId = 0; //Initially not assigned to any cluster
}

Point::Point(int id, double fit)
{
    dimensions = 1;
    pointId = id;
    double val;
    values.push_back(fit);

    clusterId = 0; //Initially not assigned to any cluster
}

Point::Point(int id, vector<double> v)
{
    dimensions = 0;
    pointId = id;
    //std::cout << "id in:" << id << '\n';
    //std::cout << "v size:" << v.size() << '\n';
    fit = v[0];
    for(int i=1; i<v.size(); i++){
       // std::cout << "v[i]:" << i << "->" << v[i] << '\n';
        values.push_back(v[i]);
        dimensions++;
        
    }
    
    
    clusterId = 0; //Initially not assigned to any cluster
    //std::cout << "*********************" << '\n';
}

int Point::getDimensions()
{
    return dimensions;
}

int Point::getCluster()
{
    return clusterId;
}

double Point::getFit()
{
    return fit;
}

int Point::getID()
{
    return pointId;
}

void Point::setCluster(int val)
{
    clusterId = val;
}

double Point::getVal(int pos)
{
    return values[pos];
}

Cluster::Cluster(int clusterId, Point centroid)
{
    this->clusterId = clusterId;
    for (int i = 0; i < centroid.getDimensions(); i++)
    {
        this->centroid.push_back(centroid.getVal(i));
    }
    this->addPoint(centroid);
}

void Cluster::addPoint(Point p)
{
    p.setCluster(this->clusterId);
    points.push_back(p);
}

bool Cluster::removePoint(int pointId)
{
    int size = points.size();

    for (int i = 0; i < size; i++)
    {
        if (points[i].getID() == pointId)
        {
            points.erase(points.begin() + i);
            return true;
        }
    }
    return false;
}

int Cluster::getId()
{
    return clusterId;
}

Point Cluster::getPoint(int pos)
{
    return points[pos];
}

int Cluster::getSize()
{
    return points.size();
}

double Cluster::getCentroidByPos(int pos)
{
    return centroid[pos];
}

void Cluster::setCentroidByPos(int pos, double val)
{
    this->centroid[pos] = val;
}

int KMeans::getNearestClusterId(Point point)
{
    double sum = 0.0, min_dist;
    int NearestClusterId;

    for (int i = 0; i < dimensions; i++)
    {
        sum += pow(clusters[0].getCentroidByPos(i) - point.getVal(i), 2.0);
    }

    min_dist = sqrt(sum);
    NearestClusterId = clusters[0].getId();

    for (int i = 1; i < K; i++)
    {
        double dist;
        sum = 0.0;

        for (int j = 0; j < dimensions; j++)
        {
            sum += pow(clusters[i].getCentroidByPos(j) - point.getVal(j), 2.0);
        }

        dist = sqrt(sum);

        if (dist < min_dist)
        {
            min_dist = dist;
            NearestClusterId = clusters[i].getId();
        }
    }

    return NearestClusterId;
}

KMeans::KMeans(int K, int iterations)
{
    this->K = K;
    this->iters = iterations;
}

void KMeans::run(double point_fits[], int size)
{
    int pointId = 0;
    vector<Point> all_points;
    for (int x = 0; x < size; x++)
    {
        Point point(pointId, point_fits[x]);
        all_points.push_back(point);
        pointId++;
    }
    total_points = all_points.size();
    dimensions = all_points[0].getDimensions();

    //Initializing Clusters
    vector<int> used_pointIds;

    for (int i = 1; i <= K; i++)
    {
        while (true)
        {
            int index = rand() % total_points;

            if (find(used_pointIds.begin(), used_pointIds.end(), index) == used_pointIds.end())
            {
                used_pointIds.push_back(index);
                all_points[index].setCluster(i);
                Cluster cluster(i, all_points[index]);
                clusters.push_back(cluster);
                break;
            }
        }
    }
    //cout<<"Clusters initialized = "<<clusters.size()<<endl<<endl;

    //cout<<"Running K-Means Clustering.."<<endl;

    int iter = 1;
    while (true)
    {
        //cout<<"Iter - "<<iter<<"/"<<iters<<endl;
        bool done = true;

        // Add all points to their nearest cluster
        for (int i = 0; i < total_points; i++)
        {
            int currentClusterId = all_points[i].getCluster();
            int nearestClusterId = getNearestClusterId(all_points[i]);

            if (currentClusterId != nearestClusterId)
            {
                if (currentClusterId != 0)
                {
                    for (int j = 0; j < K; j++)
                    {
                        if (clusters[j].getId() == currentClusterId)
                        {
                            clusters[j].removePoint(all_points[i].getID());
                        }
                    }
                }

                for (int j = 0; j < K; j++)
                {
                    if (clusters[j].getId() == nearestClusterId)
                    {
                        clusters[j].addPoint(all_points[i]);
                    }
                }
                all_points[i].setCluster(nearestClusterId);
                done = false;
            }
        }

        // Recalculating the center of each cluster
        for (int i = 0; i < K; i++)
        {
            int ClusterSize = clusters[i].getSize();

            for (int j = 0; j < dimensions; j++)
            {
                double sum = 0.0;
                if (ClusterSize > 0)
                {
                    for (int p = 0; p < ClusterSize; p++)
                        sum += clusters[i].getPoint(p).getVal(j);
                    clusters[i].setCentroidByPos(j, sum / ClusterSize);
                }
            }
        }

        if (done || iter >= iters)
        {
            //cout << "Clustering completed in iteration : " <<iter<<endl<<endl;
            break;
        }
        iter++;
    }

    int r[3];
    //Print pointIds in each cluster


    for (int i = 0; i < K; i++)
    {
        //cout<<"Points in cluster "<<clusters[i].getId()<<" : ";
        double min = clusters[i].getPoint(0).getVal(0);
        r[i] = clusters[i].getPoint(0).getID();
        for (int j = 0; j < clusters[i].getSize(); j++)
        {
            // cout<<clusters[i].getPoint(j).getVal(0)<<" ";
            //cout<<clusters[i].getPoint(j).getID()<<"\n";
            if (min > clusters[i].getPoint(j).getVal(0))
            {
                min = clusters[i].getPoint(j).getVal(0);
                r[i] = clusters[i].getPoint(j).getID();
            }
        }
    }
    this->a = r[0];
    this->b = r[1];
    this->c = r[2];
        
}


void KMeans::run(vector<double> point_fit[], int size)
{
    int pointId = 0;
    vector<Point> all_points;
    //std::cout << "xx:" << size << '\n';
    for (int x = 0; x < size; x++)
    {
        
        //std::cout << "id:" << pointId << '\n';
        Point point(pointId, point_fit[x]);
        //std::cout << "a" << '\n';
        all_points.push_back(point);
        //std::cout << "b" << '\n';
        pointId++;
        //std::cout << "d" << '\n';
        
    }
    total_points = all_points.size();
    //std::cout << "e" << '\n';
    dimensions = all_points[0].getDimensions();
    //std::cout << "f" << '\n';

    //Initializing Clusters
    vector<int> used_pointIds;

    for (int i = 1; i <= K; i++)
    {
        while (true)
        {
            int index = rand() % total_points;

            if (find(used_pointIds.begin(), used_pointIds.end(), index) == used_pointIds.end())
            {
                used_pointIds.push_back(index);
                all_points[index].setCluster(i);
                Cluster cluster(i, all_points[index]);
                clusters.push_back(cluster);
                break;
            }
        }
    }
    //cout<<"Clusters initialized = "<<clusters.size()<<endl<<endl;

    //cout<<"Running K-Means Clustering.."<<endl;
    int iter = 1;
    while (true)
    {
        //cout<<"Iter - "<<iter<<"/"<<iters<<endl;
        bool done = true;

        // Add all points to their nearest cluster
        for (int i = 0; i < total_points; i++)
        {
            int currentClusterId = all_points[i].getCluster();
            int nearestClusterId = getNearestClusterId(all_points[i]);

            if (currentClusterId != nearestClusterId)
            {
                if (currentClusterId != 0)
                {
                    for (int j = 0; j < K; j++)
                    {
                        if (clusters[j].getId() == currentClusterId)
                        {
                            clusters[j].removePoint(all_points[i].getID());
                        }
                    }
                }

                for (int j = 0; j < K; j++)
                {
                    if (clusters[j].getId() == nearestClusterId)
                    {
                        clusters[j].addPoint(all_points[i]);
                    }
                }
                all_points[i].setCluster(nearestClusterId);
                done = false;
            }
        }

        // Recalculating the center of each cluster
        for (int i = 0; i < K; i++)
        {
            int ClusterSize = clusters[i].getSize();

            for (int j = 0; j < dimensions; j++)
            {
                double sum = 0.0;
                if (ClusterSize > 0)
                {
                    for (int p = 0; p < ClusterSize; p++)
                        sum += clusters[i].getPoint(p).getVal(j);
                    clusters[i].setCentroidByPos(j, sum / ClusterSize);
                }
            }
        }

        if (done || iter >= iters)
        {
            //cout << "Clustering completed in iteration : " <<iter<<endl<<endl;
            break;
        }
        iter++;
    }
    //std::cout << "h" << '\n';
    int r = 0;
    //Print pointIds in each cluster


    int choose_cluster=0;
    std::cout << "cluster " << choose_cluster << " size: " << clusters[choose_cluster].getSize() << '\n';
    for (int i = 1; i < K; i++)
    {
        std::cout << "cluster " << i << " size: " << clusters[i].getSize() << '\n';
        if (clusters[choose_cluster].getSize() < clusters[i].getSize())
            choose_cluster = i;
    }
    
    double min = clusters[choose_cluster].getPoint(0).getFit();
    
    for (int j = 1; j < clusters[choose_cluster].getSize(); j++)
    {
        if (min > clusters[choose_cluster].getPoint(j).getFit())
        {
            min = clusters[choose_cluster].getPoint(j).getFit();
            r = clusters[choose_cluster].getPoint(j).getID();
        }
    }
    this->a = r;

     
    //std::cout << "j" << '\n';
    std::cout << "choose cluster: " << choose_cluster << '\n';

    /*
    std::cout << r << '\n';
    std::cout << min << '\n';*/
}