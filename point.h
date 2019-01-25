#pragma once

#include <math.h>
using namespace std;

class Vector{
    public:
        double i;
        double j;
        double k;
        Vector(double a, double b, double c):i(a), j(b), k(c){}
};

class Point{
    public:
        double x;
        double y;
        double z;
        Point(double a, double b, double c):x(a), y(b), z(c){}
};

class IntersectionPoint{
    Point location; // location of the intersection
    Vector normal; // normal at the location of intersection
    double rayParameter; //p0 + tp1 // tha value of t
    public:
        IntersectionPoint(Point&p, Vector& n, double r) : location(p), normal(n), rayParameter(r){}
        Point getLocation(){
            return location;
        }
        Vector getNormal(){
            return normal;
        }
        double getRayParameter(){
            return rayParameter;
        }
};

Vector getUnitVector(Vector& v){
    double norm = sqrt((v.i)*(v.i) + (v.j)*(v.j) + (v.k)*(v.k));
    Vector unitVector = v;
    unitVector.i /= norm;
    unitVector.j /= norm;
    unitVector.k /= norm;
    return unitVector;
}

vector<double> multiplyVectorsPointwise(vector<double>&v1, vector<double>&v2){
    // Both vectors should be of same size
    vector<double> result;
    for(int i=0;i<v1.size();i++){
        result.push_back(v1[i]*v2[i]);
    }
    return result;
}

vector<double> multiplyVectorDouble(const double a, const vector<double>& v){
    vector<double> result;
    for(int i=0;i<v.size();i++){
        result.push_back(v[i]*a);
    }
    return result;
}

Vector multiplyVectorDouble(const double a, const Vector& v){
    Vector result = v;
    result.i *= a;
    result.j *= a;
    result.k *= a;
    return result;
}

Point addPointVector(const Point& p, const Vector& v){
    Point pt = p;
    pt.x += v.i;
    pt.y += v.j;
    pt.z += v.k;
    return pt;
}

Vector crossProduct(Point p1, Point p2, Point p3){
    double u1 = p1.x - p2.x;
    double u2 = p1.y - p2.y;
    double u3 = p1.z - p2.z;
    double v1 = p2.x - p3.x;
    double v2 = p2.y - p3.y;
    double v3 = p2.z - p3.z;
    double uvi = u2 * v3 - v2 * u3;
    double uvj = v1 * u3 - u1 * v3;
    double uvk = u1 * v2 - v1 * u2;
    double norm = sqrt(uvi*uvi+uvj*uvj+uvk*uvk);
    Vector v(uvi/norm,uvj/norm,uvk/norm);
    return v;
}

double dotProduct(Vector p1, Vector p2){
    double uv = p1.i*p2.i + p1.j*p2.j + p1.k*p2.k;
    return uv;
}

Vector getSubtractionVector(Point p1, Point p2){
    double x = p2.x - p1.x;
    double y = p2.y - p1.y;
    double z = p2.z - p1.z;
    Vector v(x,y,z);
    return v;
} 

class Ray{
    Point source;
    Vector direction;
    public:
        Ray(Point& p, Vector& v):source(p), direction(v){}
        Point getSource(){
            return source;
        }
        Vector getDirection(){
            return direction;
        }
};