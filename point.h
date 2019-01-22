#pragma once

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

Vector MultiplyVectorDouble(const double a, const Vector& v){
    Vector result = v;
    result.i *= a;
    result.j *= a;
    result.k *= a;
    return result;
}

Point AddPointVector(const Point& p, const Vector& v){
    Point pt = p;
    pt.x += v.i;
    pt.y += v.j;
    pt.z += v.k;
    return pt;
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