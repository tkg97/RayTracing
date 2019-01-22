#pragma once

#include <bits/stdc++.h>
#include "point.h"
using namespace std;

class quadric{
    //  ax^2 + by^2 + cz^2 + 2fyz + 2gzx + 2hxy + 2px + 2qy + 2rz + d =0.
    double a, b, c, f, g, h, p, q, r, d;
    double getFunctionValue(double x1, double y1, double z1, double x2, double y2, double z2){
        return (a*x1*x2 + b*y1*y2 + c*z1*z2 + f*(y1*z2 + y2*z1) + g*(z1*x2 + z2*x1) + h*(x1*y2 + x2*y1) + 
        p*(x1+x2) + q*(y1+y2) + r*(z1+z2) + d);
    }
    public:
        quadric(double d1, double d2, double d3, double d4, double d5, double d6, double d7, double d8, double d9, double d10){
            a = d1;
            b = d2;
            c = d3;
            f = d4;
            g = d5;
            h = d6;
            p = d7;
            q = d8;
            r = d9;
            d = d10;
        }
        
        Vector getNormal(Point& pt){
            double i = a*(pt.x) + g*(pt.z) + h*(pt.y) + p; // factor of 2 is removed as normalization will be done
            double j = b*(pt.y) + f*(pt.z) + h*(pt.x) + q;
            double k = c*(pt.z) + f*(pt.y) + g*(pt.x) + r;
            double norm = sqrt(i*i + j*j + k*k);
            Vector normal(i/norm,j/norm,k/norm);
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