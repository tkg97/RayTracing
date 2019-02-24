#pragma once

#include "object.h"

class Sphere : public Object {
	double radius;
	Point center;

	double getValueA(double x1, double y1, double z1, double x2, double y2, double z2) {
		return (x2*x2 + y2 * y2 + z2 * z2);
	}

	double getValueB(double x1, double y1, double z1, double x2, double y2, double z2) {
		return 2 * (x1*x2 + y1 * y2 + z1 * z2 - (center.x)*(x2)-(center.y)*(y2)-(center.z)*(z2));
	}

	double getValueC(double x1, double y1, double z1, double x2, double y2, double z2) {
		return (x1*x1 + y1 * y1 + z1 * z1 - (center.x)*(x1 + x1) + -(center.y)*(y1 + y1) - (center.z)*(z1 + z1) +
			((center.x)*(center.x) + (center.y)*(center.y) + (center.z)*(center.z) - radius * radius));
	}

	Vector getNormal(Point p) {
		Vector normal = multiplyMatrixVector(getNormalTransformationMatrix(), getSubtractionVector(center, p));
		return getUnitVector(normal);
	}

	pair<int, int> getImageCoordinates(Point p) override {
		double theta = atan2(-(p.z - center.z), (p.x - center.x));
		double phi = acos(-(p.y - center.y) / radius);
		Material objectMat = getObjectMaterial();
		int u = (int)(((theta + M_PI) / ((2.0) * M_PI))*objectMat.getImageWidth());
		int v = (int)((phi / M_PI)*objectMat.getImageHeight());
		return { u,v };
	}

public:
	Sphere(double rad, Point pt, Material m, vector<vector<double>> t) : Object(m), radius(rad), center(pt) {
		vector<vector<double>> translationMat = formTranslationMatrix(center.x, center.y, center.z);
		setTransformations(t, translationMat);
	}

	IntersectionPoint* getIntersection(Ray r, double minThreshold) {
		r = getTransformedRay(r, getInverseVertexTransformationMatrix(), getRayTransformation());
		Point raySource = r.getSource();
		Vector rayDirection = r.getDirection();
		double a = getValueA(raySource.x, raySource.y, raySource.z, rayDirection.i, rayDirection.j, rayDirection.k);
		double b = getValueB(raySource.x, raySource.y, raySource.z, rayDirection.i, rayDirection.j, rayDirection.k);
		double c = getValueC(raySource.x, raySource.y, raySource.z, rayDirection.i, rayDirection.j, rayDirection.k);
		double t;
		if (abs(a) <= 0.00001) {
			t = -c / b;
		}
		else {
			if (b*b - 4 * a*c < 0) t = -1;
			else {
				double t1 = (-b - sqrt(b*b - 4 * a*c)) / (2 * a);
				double t2 = (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
				if (t1 >= minThreshold && t2 >= minThreshold) {
					t = min(t1, t2);
				}
				else if (t1 >= minThreshold) {
					t = t1;
				}
				else if (t2 >= minThreshold) {
					t = t2;
				}
				else t = -1;
			}
		}
		if (t < minThreshold) return nullptr;
		else {
			Point l = (addPointVector(raySource, multiplyVectorDouble(t, rayDirection)));
			Vector n = getNormal(l);
			return (new IntersectionPoint(multiplyMatrixVector(getVertexTransformation(), l), n));
		}
	}
};
