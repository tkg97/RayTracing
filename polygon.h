#pragma once
#include <bits/stdc++.h>
#include "point.h"
using namespace std;

class Polygon{
    int n;
    vector< Point > coordinates;
    public:
        Polygon(int t, vector<Point> & v){
            n=t;
            for(int i=0;i<n;i++){
               coordinates.push_back(v[i]); 
            }
        }

        // return normal to the polygon using 3 corner points
        Vector getNormal(){
             Point p1 = coordinates[0];
             Point p2 = coordinates[1];
             Point p3 = coordinates[2];
             Vector normal = crossProduct(p1,p2,p3);
             return normal;
        }
        
        // get the intersection of ray with the polygon
        Point* getIntersection(Ray r1){
             Point r0 = r1.getSource();
             Point p0 = coordinates[0];
             Vector rd = r1.getDirection();
             Vector normal = getNormal();
             double x1 = normal.i;
             double x2 = normal.j;
             double x3 = normal.k;
             double a1 = r0.x - p0.x;
             double a2 = r0.y - p0.y;
             double a3 = r0.z - p0.z;
             double b1 = rd.i;
             double b2 = rd.j;
             double b3 = rd.k;
             double t = 0-(x1*a1 + x2*a2 + x3*a3)/(x1*b1 + x2*b2 + x3*b3);
             Point intersection(r0.x+rd.i*t,r0.y+rd.j*t,r0.z+rd.k*t);
             Point *p = new Point(r0.x+rd.i*t,r0.y+rd.j*t,r0.z+rd.k*t);
             bool c = isContained(intersection);
             if(c)
                return p;   
             else
                return NULL;
        }
        
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
            double x_comp = q1.x - q2.x;
            double y_comp = q1.y - q2.y;
            double z_comp = q1.z - q2.z;
            // t1 and t2 represent the intersection point
            double t1 = (y_comp*q1.x - x_comp*q1.y-p0.x*y_comp + p0.y*x_comp)/(p1.i*y_comp - p1.j*x_comp);
            double t2 = (z_comp*q1.x - x_comp*q1.z-p0.x*z_comp + p0.z*x_comp)/(p1.i*z_comp - p1.k*x_comp);
            // return the ratio of distances of the line segments 
            double slope = (q1.x - p0.x - p1.i*t1)/(q1.x - q2.x);
            // check if t1 and t2 are same and the ratio of distances should lie between 0 and 1
            if(abs(t1-t2)<=0.0001 && slope >=-0.0001 && slope < 1.0001 ){
                return true;
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
            double norm = sqrt(u1*u1+u2*u2+u3*u3);
            Vector direction(u3/norm,u3/norm,u3/norm); 
        
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
        
            // Return true if count is odd, false otherwise 
            return count&1;  // Same as (count%2 == 1) 
        } 

        
};