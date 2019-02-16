#pragma once

#include <vector>
#include "point.h"

using namespace std;

class Viewer{
	Point eyeLocation;
    int width, height;
    double angle; // multiple of M_PI
    vector <vector <double> > transformationMatrix;
    public:
        Viewer(double f, int w, int h, vector< vector<double> >m) : angle(f), width(w), height(h), transformationMatrix(m), eyeLocation(0,0,0){
			eyeLocation = multiplyMatrixVector(transformationMatrix, eyeLocation);
		}
        Point getEyeLocation(){
			return eyeLocation;
        }
        double getAngle(){
            return angle;
        }
        int getWidth(){
            return width;
        }
        int getHeight(){
            return height;
        }
        vector<vector<double>> getTransformationMatrix(){
            return transformationMatrix;
        }

		Vector getRayDirection(double x, double y) {
			Point a(x, y, -1);
			Point b = multiplyMatrixVector(transformationMatrix, a);
			Vector dir = getSubtractionVector(eyeLocation, b);
			return dir;
		}
};