#pragma once

#include <vector>
#include "point.h"
using namespace std;

class LightSource{
    Point location;
    vector<double> intensity;
    public:
        LightSource(Point p, vector<double> d) : location(p), intensity(d){
            
        }
        Point getLocation(){
            return location;
        }
        vector<double> getIntensity(){
            return intensity;
        }
}