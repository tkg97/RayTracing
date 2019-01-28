#pragma once

#include <vector>
#include <iostream>
#include <algorithm>
#include <float.h>
#include "point.h"
using namespace std;

class Object{
    // an abstract parent class for all the objects
    vector<double> ambientCoefficient;
    vector<double> diffuseCoefficient;
    vector<double> specularCoefficient;
    double refractiveIndex;
    double phongExponent;
	double reflectionConstant;
	double refractionConstant;
    bool planarity;
    public:
        Object(vector<double> a, vector<double> d, vector<double> s, double r, double p, double c1, double c2, bool planar = false) : ambientCoefficient(a), diffuseCoefficient(d), specularCoefficient(s),
                     refractiveIndex(r), phongExponent(p), reflectionConstant(c1), refractionConstant(c2), planarity(planar){}
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
        bool isPlanar(){
            return planarity;
        }
        virtual IntersectionPoint* getIntersection(Ray r, double minThreshold) = 0;
        // minThreshold will help me handle the case for t==0
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
        Vector normal = getSubtractionVector(center, p);
        return getUnitVector(normal);
    }

    public:
        Sphere(double rad, Point pt, vector<double> a, vector<double> d, vector<double> s, double r, double p, double c1, double c2): Object(a, d, s, r, p, c1, c2), radius(rad), center(pt){}

        IntersectionPoint* getIntersection(Ray r, double minThreshold){
            // cout << "i am called" << endl;
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
    vector<Vector> normals;
    vector<Point> referencePoints;

    void calculateAllNormals(){
        normals.push_back(crossProduct(coordinates[2], coordinates[1], coordinates[0]));
        normals.push_back(crossProduct(coordinates[0], coordinates[4], coordinates[7]));
        normals.push_back(crossProduct(coordinates[4], coordinates[5], coordinates[6]));
        normals.push_back(crossProduct(coordinates[5], coordinates[1], coordinates[2]));
        normals.push_back(crossProduct(coordinates[5], coordinates[1], coordinates[0]));
        normals.push_back(crossProduct(coordinates[2], coordinates[3], coordinates[7]));
    }

    void storeReferencePoints(){
        referencePoints.push_back(coordinates[2]);
        referencePoints.push_back(coordinates[0]);
        referencePoints.push_back(coordinates[4]);
        referencePoints.push_back(coordinates[5]);
        referencePoints.push_back(coordinates[5]);
        referencePoints.push_back(coordinates[2]);
    }
    
    public:
        Box(vector<Point> v, vector<double> a, vector<double> d, vector<double> s, double r, double p, double c1, double c2) : Object(a, d, s, r, p, c1, c2), coordinates(v){
            calculateAllNormals();
            storeReferencePoints();
        }

        IntersectionPoint* getIntersection(Ray r, double minThreshold){
            Point raySource = r.getSource();
            Vector rayDirection = r.getDirection(); 
            double tEnterMax = DBL_MIN, tLeaveMin = DBL_MAX;
            Vector nEnter(0, 0, 0); // Initialised to be a zero vector // normal
			Vector nLeave(0, 0, 0); // Initialised to be a zero vector // normal
            for(int i=0;i<6;i++){
                if(dotProduct(rayDirection, normals[i]) == 0){
                // if ray is parallel to certain face and it is not contained in the face, then its probably useless
                // And the case where it is contained is simply avoided
                    continue;
                }
                if(dotProduct(rayDirection, normals[i]) <0){
                    //ray is entering
                    double t = -dotProduct(getSubtractionVector(referencePoints[i], raySource), normals[i])/(dotProduct(rayDirection, normals[i]));
                    if(t>tEnterMax){
                        tEnterMax = t;
                        nEnter = normals[i];
                    }
                }
                else{
                    //ray is leaving
                    double t = -dotProduct(getSubtractionVector(referencePoints[i], raySource), normals[i])/(dotProduct(rayDirection, normals[i]));
					if (t < tLeaveMin) {
						tLeaveMin = t;
						nLeave = normals[i];
					}
                }
            }
            if(tEnterMax > tLeaveMin) return nullptr;
            else{
				if (tEnterMax >= minThreshold) {
					Point l = (addPointVector(raySource, multiplyVectorDouble(tEnterMax, rayDirection)));
					return (new IntersectionPoint(l, nEnter, tEnterMax));
				}
				else {
					Point l = (addPointVector(raySource, multiplyVectorDouble(tLeaveMin, rayDirection)));
					return (new IntersectionPoint(l, nLeave, tLeaveMin));
				}
            }
        }
};

class quadric : public Object{
    //  ax^2 + by^2 + cz^2 + 2fyz + 2gzx + 2hxy + 2px + 2qy + 2rz + d = 0.
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
        return normal;
    }

    public:
        quadric(double d1, double d2, double d3, double d4, double d5, double d6, double d7, double d8, double d9, double d10,
			vector<double> amb, vector<double> d, vector<double> s, double r, double p, double c1, double c2, bool planar = false) : Object(amb, d, s, r, p, c1, c2, planar),
            a(d1), b(d2), c(d3), f(d4), g(d5), h(d6), p(d7), q(d8), r(d9), d(d10){}

        IntersectionPoint* getIntersection(Ray r, double minThreshold){
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

class Polygon : public Object{
    int n; //vertexCount
    vector< Point > coordinates;
    //returns true if q lies on line segment p-r
    bool onSegment(Point p, Point q, Point r) 
    { 
        if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) && 
            q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y) && 
            q.z <= max(p.z, r.z) && q.z >= min(p.z, r.z) ) 
            return true; 
        return false; 
    } 

    // To find orientation of ordered triplet (p, q, r). 
    // The function returns following values 
    // 0 --> p, q and r are colinear 
    // 1 --> if not
    int orientation(Point p, Point q, Point r) 
    { 
        double val1 = (q.y - p.y) * (r.x - q.x) - 
                (q.x - p.x) * (r.y - q.y); 
        double val2 = (q.z - p.z) * (r.x - q.x) - 
                (q.x - p.x) * (r.z - q.z);  
        if (abs(val1) <= 0.0001 && abs(val2)<=0.0001 ) return 0;  // colinear 
        // return (val > 0)? 1: 2; // clock or counterclock wise 
        return  1;
    } 

    // The function that returns true if the ray 'p0 + p1*t' 
    // and the line segment 'q1q2' intersect. 
    bool doIntersect(Point p0, Point q1, Point q2, Vector p1) 
    { 
        double a1 = p1.i , a2= p1.j , a3 = p1.k;
        double b1 = q1.x - q2.x;
        double b2 = q1.y - q2.y;
        double b3 = q1.z - q2.z;
        double c1 = q1.x - p0.x , c2 = q1.y - p0.y , c3  = q1.z - p0.z;
        // t1 and t2 represent the intersection point
        // cout<<x_comp<<" "<<y_comp<<" "<<z_comp<<" "<<p1.i<<" "<<p1.j<<" "<<p1.k<<endl;
        // cout<<p1.i*y_comp - p1.j*x_comp<<" "<<p1.i*z_comp - p1.k*x_comp<<" "<<p1.j*z_comp - p1.k*y_comp<<" ";
        double t1,t2,t3;
        double r1,r2,r3;
        // if(abs(t1-t2)>0.0001 || abs(t2-t3)>0.0001 || abs(t1-t3)>0.0001) return false;
        if(abs(a1*b2 - a2*b1) >= 0.0001){
            t1 = (c1*b2 - b1*c2)/(a1*b2 - a2*b1);
            r1 = (-c1*a2 + a1*c2)/(a1*b2 - a2*b1);
            if(abs(p0.z + p1.k*t1 - (q1.z + r1*(q2.z - q1.z)))<=0.001 && r1 <= 1.0001 && r1 >=-0.0001 && t1>=-0.0001)
               return true;
        }
        else if(abs(a2*b3 - a3*b2) >= 0.0001){
            t1 = (c2*b3 - b2*c3)/(a2*b3 - a3*b2);
            r1 = (-c2*a3 + a2*c3)/(a2*b3 - a3*b2);
            if(abs(p0.x + p1.i*t1 - (q1.x + r1*(q2.x - q1.x)))<=0.001 && r1 <= 1.0001 && r1 >=-0.0001 && t1>=-0.0001)
               return true;
        }
        else if(abs(a1*b3 - a3*b1) >= 0.0001){
            t1 = (c1*b3 - b1*c3)/(a1*b3 - a2*b3);
            r1 = (-c1*a3 + a1*c3)/(a1*b3 - a3*b1);
            if(abs(p0.y + p1.j*t1 - (q1.y + r1*(q2.y - q1.y)))<=0.001 && r1 <= 1.0001 && r1 >=-0.0001 && t1>=-0.0001)
               return true;
        }
        else{
            bool f = onSegment( q1 ,p0 , q2);
            if(f==true) return true;
        }
        return false;
    } 

    // Returns true if the point p lies inside the polygon[] with n vertices 
    bool isContained(Point p) 
    { 
        // There must be at least 3 vertices in polygon[] 
        if (n < 3)  return false; 
        
        // Create a direction for a ray from p parallel to some edge
        Point p1 = coordinates[0];
        Point p2 = coordinates[1];
        double u1 = p1.x - p2.x;
        double u2 = p1.y - p2.y;
        double u3 = p1.z - p2.z;
        // cout<<u1<<" "<<u2<<" "<<u3<<endl;
        double norm = sqrt(u1*u1+u2*u2+u3*u3);
        Vector direction(u1/norm,u2/norm,u3/norm); 
        
        // Count intersections of the above line with sides of polygon 
        int count = 0, i = 0; 
        do
        { 
            int next = (i+1)%n; 
    
            // Check if the line segment from 'p' in direction given by 'direction' intersects 
            // with the line segment from 'coordinates[i]' to 'coordinates[next]' 
            if (doIntersect(p,coordinates[i], coordinates[next], direction)) 
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
        // cout<<count<<" ";
        // Return true if count is odd, false otherwise 
        return count&1;  // Same as (count%2 == 1) 
    }

    // return normal to the polygon using 3 corner points
    Vector getNormal(){
        Point p1 = coordinates[0];
        Point p2 = coordinates[1];
        Point p3 = coordinates[2];
        Vector normal = crossProduct(p1,p2,p3);
        return normal;
    }

    public:
        Polygon(int t, vector<Point> v, vector<double> a, vector<double> d, vector<double> s, double r, double p, double c1, double c2) : Object(a, d, s, r, p, c1, c2, true){
            n=t;
            for(int i=0;i<n;i++){
               coordinates.push_back(v[i]); 
            }
        }
        
        // get the intersection of ray with the polygon
        IntersectionPoint* getIntersection(Ray r1, double minThreshold){
             Point r0 = r1.getSource();
             Point p0 = coordinates[0];
             Vector rd = r1.getDirection();
             Vector normal = getNormal();
             double a1 = r0.x - p0.x;
             double a2 = r0.y - p0.y;
             double a3 = r0.z - p0.z;
             Vector v(a1,a2,a3);
             if(abs(dotProduct(normal,rd))<=0.00001){
                  return nullptr;
             }
             double t = -(dotProduct(normal,v))/(dotProduct(normal,rd));
			 if (t < minThreshold) return nullptr;
             
			 Point intersection = addPointVector(r0, multiplyVectorDouble(t,rd));          
             if(isContained(intersection)){
                return (new IntersectionPoint(intersection, getNormal(), t));  
             } 
             else
                return nullptr;
        }     
};