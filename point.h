#pragma once

#include <cmath>
#include <vector>
using namespace std;

class Vector{
    public:
        double i;
        double j;
        double k;
        Vector(double a, double b, double c):i(a), j(b), k(c){}
        Vector operator +(const Vector&v){
            Vector result = v;
            result.i += this->i;
            result.j += this->j;
            result.k += this->k;
            return result;
        }
};

class Point{
    public:
        double x;
        double y;
        double z;
        Point(double a, double b, double c):x(a), y(b), z(c){}
};

Vector getUnitVector(Vector v) {
	double norm = sqrt((v.i)*(v.i) + (v.j)*(v.j) + (v.k)*(v.k));
	v.i /= norm;
	v.j /= norm;
	v.k /= norm;
	return v;
}

Vector multiplyVectorDouble(double a, Vector v) {
	v.i *= a;
	v.j *= a;
	v.k *= a;
	return v;
}

vector<double> multiplyVectorsPointwise(vector<double>v1, vector<double>v2) {
	// Both vectors should be of same size
	for (int i = 0;i < v1.size();i++) {
		v1[i] *= v2[i];
	}
	return v1;
}

vector<double> multiplyVectorDouble(double a, vector<double> v) {
	for (int i = 0;i < v.size();i++) {
		v[i] *= a;
	}
	return v;
}

Point addPointVector(Point p, Vector v) {
	p.x += v.i;
	p.y += v.j;
	p.z += v.k;
	return p;
}

Vector crossProduct(Point p1, Point p2, Point p3) {
	// forms two vectors u(p2p1) and v(p3p2) and returns u cross v
	double u1 = p1.x - p2.x;
	double u2 = p1.y - p2.y;
	double u3 = p1.z - p2.z;
	double v1 = p2.x - p3.x;
	double v2 = p2.y - p3.y;
	double v3 = p2.z - p3.z;
	double uvi = u2 * v3 - v2 * u3;
	double uvj = v1 * u3 - u1 * v3;
	double uvk = u1 * v2 - v1 * u2;
	double norm = sqrt(uvi*uvi + uvj * uvj + uvk * uvk);
	Vector v(uvi / norm, uvj / norm, uvk / norm);
	return v;
}

double dotProduct(Vector p1, Vector p2) {
	double uv = p1.i*p2.i + p1.j*p2.j + p1.k*p2.k;
	return uv;
}

Vector getSubtractionVector(Point p1, Point p2) {
	// vector p1p2 (p2 - p1)
	double x = p2.x - p1.x;
	double y = p2.y - p1.y;
	double z = p2.z - p1.z;
	Vector v(x, y, z);
	return v;
}

double getNorm(Vector v) {
	return sqrt(v.i*v.i + v.j*v.j + v.k*v.k);
}

class Ray{
    Point source;
    Vector direction;
    public:
        Ray(Point p, Vector v) : source(p), direction(getUnitVector(v)){}
        Point getSource(){
            return source;
        }
        Vector getDirection(){
            return direction;
        }
};

class IntersectionPoint{
    Point location; // location of the intersection
    Vector normal; // normal at the location of intersection
    public:
        IntersectionPoint(Point p, Vector n) : location(p), normal(n){}
        Point getLocation(){
            return location;
        }
        Vector getNormal(){
            return normal;
        }
};

struct RayData{
	int screenI, screenJ;
	Point start, end;

	RayData(int a, int b, Point p, Point q) : screenI(a), screenJ(b), start(p), end(q) {}
};
