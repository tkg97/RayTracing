#pragma once

#include <vector>
#include "point.h"

using namespace std;

class Viewer{
    Point eyeLocation;
    int width, height;
    double angle;
    public:
        Viewer(Point& e, double f, int w, int h) : eyeLocation(e),angle(f), width(w), height(h){}
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
};