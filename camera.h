#pragma once

#include <vector>
#include "point.h"

using namespace std;

class Viewer{
    Point eyeLocation;
    int width, height;
    double angle;
    vector<Vector> orthoDirections ; // defining axis system wrt to image plane // third one is the viewingDirection. 
    public:
        Viewer(Point& e, double f, int w, int h, vector<Vector>& v) : eyeLocation(e),angle(f), width(w), height(h), orthoDirections(v){
            // normalizing orthoDirections for future use
            for(int d=0;d<3;d++){
                double norm = sqrt((orthoDirections[d].i * orthoDirections[d].i) 
                                + (orthoDirections[d].j * orthoDirections[d].j) + 
                                (orthoDirections[d].k * orthoDirections[d].k));
                orthoDirections[d].i /= norm;
                orthoDirections[d].j /= norm;
                orthoDirections[d].k /= norm;
            }
        }
        Point getEyeLocation(){
            return eyeLocation;
        }
        double getAngle(){
            return angle;
        }
        int getWidth(){
            return width;
        }
        int getHeight(){
            return height;
        }
        vector<Vector> getOrthoDirections(){
            return orthoDirections;
        }
};