#pragma once 

#include<vector>
#include "scene.h"
#include "point.h"
#define  M_PI 3.14
using namespace std;

// multiply a matrix(4*4) with an augmented point to get the point in the new coordinate system
Point multiplyMatrix(vector<vector<double>> cameraToWorld,Point p){
    double u1 = p.x*cameraToWorld[0][0] + p.y*cameraToWorld[1][0] + p.z*cameraToWorld[2][0] + cameraToWorld[3][0];
    double u2 = p.x*cameraToWorld[0][1] + p.y*cameraToWorld[1][1] + p.z*cameraToWorld[2][1] + cameraToWorld[3][1];
    double u3 = p.x*cameraToWorld[0][2] + p.y*cameraToWorld[1][2] + p.z*cameraToWorld[2][2] + cameraToWorld[3][2];
    Point r(u1,u2,u3);
    return r;
}
int main(){
    vector<vector<double>> cameraToWorld; 
    double rgb[width][height]; 
    double scale = tan(angle * M_PI / 180); 
    double imageAspectRatio = width / height; 
    Point o(0,0,0);
    Point origin = multiplyMatrix(cameraToWorld, o);
    Scene s;
    for (int j = 0; j < height; ++j) { 
        for (int i = 0; i < width; ++i) { 
            float x = (2 * ((i + 0.5) / width) - 1) * imageAspectRatio * scale; 
            float y = (1 - 2 * ((j + 0.5) / height)) * scale; 
            Point a(x,y,-1);
            Point b = multiplyMatrix(cameraToWorld ,a );
            Vector dir = getSubtractionVector(origin,b); 
            Ray r(origin , dir);
            vector<double> rgb = s.getIlumination(r, -1); 
        } 
    } 


    return 0;
}


