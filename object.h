#pragma once

#include <vector>
#include <iostream>
#include <algorithm>
#include <float.h>
#include "point.h"
using namespace std;

class Material{
        vector<double> ambientCoefficient;
        vector<double> diffuseCoefficient;
        vector<double> specularCoefficient;
        double refractiveIndex;
        double phongExponent;
        double reflectionConstant;
        double refractionConstant;
    public:
        Material(vector<double> a, vector<double> d, vector<double> s, double r, double p, double c1, double c2) : ambientCoefficient(a), diffuseCoefficient(d), specularCoefficient(s),
                        refractiveIndex(r), phongExponent(p), reflectionConstant(c1), refractionConstant(c2){}
        vector<double> getAmbientCoefficient(){
            return ambientCoefficient;
        }
        vector<double> getDiffusionCoefficeint(){
            return diffuseCoefficient;
        }
        vector<double> getSpecularCoefficient(){
            return specularCoefficient;
        }
        double getRefractiveIndex(){
            return refractiveIndex;
        }
        double getPhongExponent(){
            return phongExponent;
        }
        double getReflectionConstant(){
            return reflectionConstant;
        }
        double getRefractionConstant(){
            return refractionConstant;
        }
};

class Object{
    // an abstract parent class for all the objects
        Material objectMaterial;
        bool planarity;
		vector<vector<double>> vertexTransformation;
		vector<vector<double>> normalTransformation;
		vector<vector<double>> rayTransformation;
    public:
        Object(Material m, vector<vector<double>> t, bool planar = false) : objectMaterial(m), planarity(planar), vertexTransformation(t)
		{
			rayTransformation = getMatrixInverse(vertexTransformation);
			normalTransformation = getMatrixTranspose(rayTransformation);
		}
        vector<double> getAmbientCoefficient(){
            return objectMaterial.getAmbientCoefficient();
        }
        vector<double> getDiffusionCoefficeint(){
            return objectMaterial.getDiffusionCoefficeint();
        }
        vector<double> getSpecularCoefficient(){
            return objectMaterial.getSpecularCoefficient();
        }
        double getRefractiveIndex(){
            return objectMaterial.getRefractiveIndex();
        }
        double getPhongExponent(){
            return objectMaterial.getPhongExponent();
        }
        double getReflectionConstant(){
            return objectMaterial.getReflectionConstant();
        }
        double getRefractionConstant(){
            return objectMaterial.getRefractionConstant();
        }
        bool isPlanar(){
            return planarity;
        }
		vector<vector<double>> getVertexTransformationMatrix() {
			return vertexTransformation;
		}
		vector<vector<double>> getNormalTransformationMatrix() {
			return normalTransformation;
		}
		vector<vector<double>> getRayTransformation() {
			return rayTransformation;
		}
        virtual IntersectionPoint* getIntersection(Ray r, double minThreshold) = 0;
        // minThreshold will help me handle the case for t==0 avoidance
};

class Polygon : public Object {
	int vertexCount; //vertexCount
	vector< Point > coordinates;
	Vector normal;
	//returns true if q lies on line segment p-r
	bool onSegment(Point p, Point q, Point r)
	// first collinearity has to be checked
	{
		if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
			q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y) &&
			q.z <= max(p.z, r.z) && q.z >= min(p.z, r.z))
			return true;
		return false;
	}

	// To find orientation of ordered triplet (p, q, r). 
	// The function returns following values 
	// 0 --> p, q and r are colinear 
	// 1 --> if not
	int orientation(Point p, Point q, Point r)
	{
		Vector crossProductVal = crossProduct(p, q, r);
		if (abs(crossProductVal.i) <= 0.0001 && abs(crossProductVal.j) <= 0.0001 && abs(crossProductVal.k) <= 0.0001) return 0;  // collinear  
		return  1;
	}

	// The function that returns true if the ray 'p0 + p1*t' 
	// and the line segment 'q1q2' intersect. 
	bool doIntersect(Point p0, Point q1, Point q2, Vector p1)
	{
		double a1 = p1.i, a2 = p1.j, a3 = p1.k;
		double b1 = q1.x - q2.x;
		double b2 = q1.y - q2.y;
		double b3 = q1.z - q2.z;
		double c1 = q1.x - p0.x, c2 = q1.y - p0.y, c3 = q1.z - p0.z;
		// t1 and t2 represent the intersection point
		
		double t1, t2, t3;
		double r1, r2, r3;
		
		if (abs(a1*b2 - a2 * b1) >= 0.0001) {
			t1 = (c1*b2 - b1 * c2) / (a1*b2 - a2 * b1);
			r1 = (-c1 * a2 + a1 * c2) / (a1*b2 - a2 * b1);
			if (abs(p0.z + p1.k*t1 - (q1.z + r1 * (q2.z - q1.z))) <= 0.001 && r1 <= 1.0001 && r1 >= -0.0001 && t1 >= -0.0001)
				return true;
		}
		else if (abs(a2*b3 - a3 * b2) >= 0.0001) {
			t1 = (c2*b3 - b2 * c3) / (a2*b3 - a3 * b2);
			r1 = (-c2 * a3 + a2 * c3) / (a2*b3 - a3 * b2);
			if (abs(p0.x + p1.i*t1 - (q1.x + r1 * (q2.x - q1.x))) <= 0.001 && r1 <= 1.0001 && r1 >= -0.0001 && t1 >= -0.0001)
				return true;
		}
		else if (abs(a1*b3 - a3 * b1) >= 0.0001) {
			t1 = (c1*b3 - b1 * c3) / (a1*b3 - a3 * b1);
			r1 = (-c1 * a3 + a1 * c3) / (a1*b3 - a3 * b1);
			if (abs(p0.y + p1.j*t1 - (q1.y + r1 * (q2.y - q1.y))) <= 0.001 && r1 <= 1.0001 && r1 >= -0.0001 && t1 >= -0.0001)
				return true;
		}
		else {
			bool c = !(orientation(q1, p0, q2));
			bool f = onSegment(q1, p0, q2);
			return c & f;
		}
		return false;
	}

	// Returns true if the point p lies inside the polygon[] with vertexCount vertices 
	bool isContained(Point p)
	{
		// There must be at least 3 vertices in polygon[] 
		if (vertexCount < 3)  return false;

		// Create a direction for a ray from p parallel to some edge
		Point p1 = coordinates[0];
		Point p2 = coordinates[1];
		double u1 = p1.x - p2.x;
		double u2 = p1.y - p2.y;
		double u3 = p1.z - p2.z;
		double norm = sqrt(u1*u1 + u2 * u2 + u3 * u3);
		Vector direction(u1 / norm, u2 / norm, u3 / norm);

		// Count intersections of the above line with sides of polygon 
		int count = 0, i = 0;
		do
		{
			int next = (i + 1) % vertexCount;

			// Check if the line segment from 'p' in direction given by 'direction' intersects 
			// with the line segment from 'coordinates[i]' to 'coordinates[next]' 
			if (doIntersect(p, coordinates[i], coordinates[next], direction))
			{
				// If the point 'p' is colinear with line segment 'i-next', 
				// then check if it lies on segment. If it lies, return true, 
				// otherwise false 
				if (orientation(coordinates[i], p, coordinates[next]) == 0)
					return onSegment(coordinates[i], p, coordinates[next]);

				count++;
			}
			i = next;
		} while (i != 0);
		// Return true if count is odd, false otherwise 
		return count & 1;  // Same as (count%2 == 1) 
	}

	// return normal to the polygon using 3 corner points
	Vector getNormal() {
		return getUnitVector(multiplyMatrix(getNormalTransformationMatrix(), normal));
	}

public:
	Polygon(int n, vector<Point> v, Material m, vector<vector<double>> t) : Object(m,t,true), vertexCount(n), coordinates(v), normal(0, 0, 0) {
		Point p1 = coordinates[0];
		Point p2 = coordinates[1];
		Point p3 = coordinates[2];
		normal = crossProduct(p1, p2, p3);
	}

	// get the intersection of ray with the polygon
	IntersectionPoint* getIntersection(Ray r1, double minThreshold) {
		r1 = getTransformedRay(r1, getRayTransformation());
		Point r0 = r1.getSource();
		Point p0 = coordinates[0];
		Vector rd = r1.getDirection();
		Vector normal = getNormal();
		double a1 = r0.x - p0.x;
		double a2 = r0.y - p0.y;
		double a3 = r0.z - p0.z;
		Vector v(a1, a2, a3);
		if (abs(dotProduct(normal, rd)) <= 0.00001) {
			return nullptr;
		}
		double t = -(dotProduct(normal, v)) / (dotProduct(normal, rd));
		if (t < minThreshold) return nullptr;

		Point intersection = addPointVector(r0, multiplyVectorDouble(t, rd));
		if (isContained(intersection)) {
			return (new IntersectionPoint(intersection, getNormal(), t));
		}
		else
			return nullptr;
	}
};

class Sphere : public Object{
    double radius;
    Point center;

    double getValueA(double x1, double y1, double z1, double x2, double y2, double z2){
        return (x2*x2 + y2*y2 + z2*z2);
    }

    double getValueB(double x1, double y1, double z1, double x2, double y2, double z2){
        return 2*(x1*x2 + y1*y2 + z1*z2 - (center.x)*(x2) - (center.y)*(y2) - (center.z)*(z2));
    }

    double getValueC(double x1, double y1, double z1, double x2, double y2, double z2){
        return (x1*x1 + y1*y1 + z1*z1 - (center.x)*(x1+x1) + - (center.y)*(y1+y1) - (center.z)*(z1+z1) + 
                ((center.x)*(center.x) + (center.y)*(center.y) + (center.z)*(center.z) - radius*radius));
    }

    Vector getNormal(Point p){
        Vector normal = multiplyMatrix(getNormalTransformationMatrix() ,getSubtractionVector(center, p));
        return getUnitVector(normal);
    }

    public:
        Sphere(double rad, Point pt, Material m, vector<vector<double>> t): Object(m,t), radius(rad), center(pt){}

        IntersectionPoint* getIntersection(Ray r, double minThreshold){
			r = getTransformedRay(r, getRayTransformation());
            Point raySource = r.getSource();
            Vector rayDirection = r.getDirection();
            double a = getValueA(raySource.x, raySource.y, raySource.z, rayDirection.i, rayDirection.j, rayDirection.k);
            double b = getValueB(raySource.x, raySource.y, raySource.z, rayDirection.i, rayDirection.j, rayDirection.k);
            double c = getValueC(raySource.x, raySource.y, raySource.z, rayDirection.i, rayDirection.j, rayDirection.k);
            double t;
            if(abs(a)<=0.00001){
                t = -c/b;
            }
            else{
                if(b*b - 4*a*c<0) t = -1;
                else{
                    double t1 = (-b - sqrt(b*b - 4*a*c))/(2*a);
                    double t2 = (-b + sqrt(b*b - 4*a*c))/(2*a);
                    if(t1>=minThreshold && t2>=minThreshold){
                        t = min(t1, t2);
                    }
                    else if(t1>=minThreshold){
                        t = t1;
                    }
                    else if(t2>=minThreshold){
                        t = t2;
                    }
                    else t = -1;
                }
            }
            if(t<minThreshold) return nullptr;
            else{
                Point l = (addPointVector(raySource, multiplyVectorDouble(t, rayDirection)));
                Vector n = getNormal(l);
                return (new IntersectionPoint(l, n, t));
            }
        }
};

class Box : public Object{
    // Assumes a rectilinear box, first four in anticlockwise order, then other four
    /*
        2 1 ---- 6 5
        | |      | |
        3 4 ---- 7 8

           5
        1 2 3 4
           6
        We will save normals for each surface in this notation
    */

    vector<Point> coordinates;
	vector<Polygon> polygonFaces;

	void storeAllPolygons(const Material &m, const vector<vector<double>> &t) {
		Polygon p1(4, { coordinates[0], coordinates[1], coordinates[2], coordinates[3] }, m, t);
		Polygon p2(4, {coordinates[7], coordinates[4], coordinates[0], coordinates[3]}, m, t);
		Polygon p3(4, {coordinates[6], coordinates[5], coordinates[4], coordinates[7]}, m, t);
		Polygon p4(4, { coordinates[2], coordinates[1], coordinates[5], coordinates[6] }, m, t);
		Polygon p5(4, { coordinates[5], coordinates[1], coordinates[0], coordinates[4] }, m, t);
		Polygon p6(4, { coordinates[7], coordinates[3], coordinates[2], coordinates[6] }, m, t);
		polygonFaces = vector<Polygon>({ p1, p2, p3, p4, p5, p6 });
	}
    
    public:
        Box(vector<Point> v, Material m, vector<vector<double> > t) : Object(m, t), coordinates(v){
			storeAllPolygons(m, t);
        }

		IntersectionPoint* getIntersection(Ray r, double minThreshold) {
			IntersectionPoint* minIntersectionPoint = nullptr;
			for (int i = 0;i < 6;i++) {
				IntersectionPoint* intersection = polygonFaces[i].getIntersection(r, minThreshold);
				if (intersection != nullptr) {
					if (minIntersectionPoint == nullptr) minIntersectionPoint = intersection;
					else if (minIntersectionPoint->getRayParameter() > intersection->getRayParameter()) minIntersectionPoint = intersection;
				}
			}
			return minIntersectionPoint;
		}
};

class Quadric : public Object{
    //  ax^2 + by^2 + cz^2 + 2fyz + 2gzx + 2hxy + 2px + 2qy + 2rz + d = 0.
    // Assumed to be hollow from inside
    double a, b, c, f, g, h, p, q, r, d;

    double getValueA(double x1, double y1, double z1, double x2, double y2, double z2){
        return (a*x2*x2 + b*y2*y2 + c*z2*z2 + f*(y2*z2 + y2*z2) + g*(z2*x2 + z2*x2) + h*(x2*y2 + x2*y2));
    }

    double getValueB(double x1, double y1, double z1, double x2, double y2, double z2){
        return 2*(a*x1*x2 + b*y1*y2 + c*z1*z2 + f*(y1*z2 + y2*z1) + g*(z1*x2 + z2*x1) + h*(x1*y2 + x2*y1) + p*(x2) + q*(y2) + r*(z2));
    }

    double getValueC(double x1, double y1, double z1, double x2, double y2, double z2){
        return (a*x1*x1 + b*y1*y1 + c*z1*z1 + f*(y1*z1 + y1*z1) + g*(z1*x1 + z1*x1) + h*(x1*y1 + x1*y1) + 
                p*(x1+x1) + q*(y1+y1) + r*(z1+z1) + d);
    }

    Vector getNormal(Point pt){
        double i = a*(pt.x) + g*(pt.z) + h*(pt.y) + p; // factor of 2 is removed as normalization will be done
        double j = b*(pt.y) + f*(pt.z) + h*(pt.x) + q;
        double k = c*(pt.z) + f*(pt.y) + g*(pt.x) + r;
        double norm = sqrt(i*i + j*j + k*k);
        Vector normal(i/norm,j/norm,k/norm);
        return getUnitVector(multiplyMatrix(getNormalTransformationMatrix(),normal));
    }

    public:
        Quadric(double d1, double d2, double d3, double d4, double d5, double d6, 
			double d7, double d8, double d9, double d10, Material m, vector<vector<double> > t)
        : Object(m,t),a(d1), b(d2), c(d3), f(d4), g(d5), h(d6), p(d7), q(d8), r(d9), d(d10){}

        IntersectionPoint* getIntersection(Ray r, double minThreshold){
			r = getTransformedRay(r, getRayTransformation());
            Point raySource = r.getSource();
            Vector rayDirection = r.getDirection();
            double a = getValueA(raySource.x, raySource.y, raySource.z, rayDirection.i, rayDirection.j, rayDirection.k);
            double b = getValueB(raySource.x, raySource.y, raySource.z, rayDirection.i, rayDirection.j, rayDirection.k);
            double c = getValueC(raySource.x, raySource.y, raySource.z, rayDirection.i, rayDirection.j, rayDirection.k);
            double t;
            if(abs(a)<=0.00001){
                t = -c/b;
            }
            else{
                if(b*b - 4*a*c<0) t = -1;
                else{
                    double t1 = (-b - sqrt(b*b - 4*a*c))/(2*a);
                    double t2 = (-b + sqrt(b*b - 4*a*c))/(2*a);
                    if(t1>=minThreshold && t2>=minThreshold){
                        t = min(t1, t2);
                    }
                    else if(t1>=minThreshold){
                        t = t1;
                    }
                    else if(t2>=minThreshold){
                        t = t2;
                    }
                    else t = -1;
                }
            }
            if(t<minThreshold) return nullptr;
            else{
                Point l = (addPointVector(raySource, multiplyVectorDouble(t, rayDirection)));
                Vector n = getNormal(l);
                return (new IntersectionPoint(l, n, t));
            }
        }
};