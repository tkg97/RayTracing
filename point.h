#pragma once
class Point{
    public:
        double x;
        double y;
        double z;
        Point(double a, double b, double c):x(a), y(b), z(c){}
};

class Vector{
    public:
        double i;
        double j;
        double k;
        Vector(double a, double b, double c):i(a), j(b), k(c){}
};

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