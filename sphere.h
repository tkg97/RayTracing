#pragma once

#include <bits/stdc++.h>
#include "point.h"
using namespace std;

class Sphere{
    double radius;
    Point center;

    double getFunctionValue(double x1, double y1, double z1, double x2, double y2, double z2){
        return (x1*x2 + y1*y2 + z1*z2 - (center.x)*(x1+x2) - (center.y)*(y1+y2) - (center.z)*(z1+z2) + 
        ((center.x)*(center.x) + (center.y)*(center.y) + (center.z)*(center.z) - radius*radius));
    }
    public:
        Sphere(double r, Point& p):radius(r), center(p){}
        
        Vector getNormal(Point p){
            double i = p.x - center.x;
            double j = p.y - center.y;
            double k = p.z - center.z;
            double norm = sqrt(i*i + j*j + k*k);
            Vector normal(i/norm, j/norm, k/norm);
            return normal;
        }

        Point* getIntersection(Ray r){
            Point raySource = r.getSource();
            Vector rayDirection = r.getDirection();
            double a = getFunctionValue(rayDirection.i, rayDirection.j, rayDirection.k, rayDirection.i, rayDirection.j, rayDirection.k);
            double b = 2*getFunctionValue(rayDirection.i, rayDirection.j, rayDirection.k, raySource.x, raySource.y, raySource.z);
            double c = getFunctionValue(raySource.x, raySource.y, raySource.z, raySource.x, raySource.y, raySource.z);
            double t;
            if(a==0){
                t = -c/b;
            }
            else{
                if(b*b - 4*a*c<0) t = -1;
                else{
                    double t1 = (-b - sqrt(b*b - 4*a*c))/(2*a);
                    double t2 = (-b + sqrt(b*b - 4*a*c))/(2*a);
                    if(t1>=0 && t2>=0){
                        t = min(t1, t2);
                    }
                    else if(t1>=0){
                        t = t1;
                    }
                    else if(t2>=0){
                        t = t2;
                    }
                    else t = -1;
                }
            }
            if(t<0) return nullptr;
            else{
                return &(AddPointVector(raySource, MultiplyVectorDouble(t, rayDirection)));
            }
        }
};