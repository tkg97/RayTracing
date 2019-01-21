#pragma once

#include <bits/stdc++.h>
#include "point.h"
using namespace std;

class Sphere{
    double radius;
    Point center;
    public:
        Sphere(double r, Point& p):radius(r), center(p){}
        Vector getNormal(Point p){
            Vector normal(p.x, p.y, p.z);
            return normal;
        }
};