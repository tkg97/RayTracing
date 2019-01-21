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