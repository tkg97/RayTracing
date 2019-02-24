#pragma once

#include "object.h"

class Quadric : public Object {
	//  ax^2 + by^2 + cz^2 + 2fyz + 2gzx + 2hxy + 2px + 2qy + 2rz + d = 0.
	// Assumed to be hollow from inside
	double a, b, c, f, g, h, p, q, r, d;
	Point center;

	double getValueA(double x1, double y1, double z1, double x2, double y2, double z2) {
		return (a*x2*x2 + b * y2*y2 + c * z2*z2 + f * (y2*z2 + y2 * z2) + g * (z2*x2 + z2 * x2) + h * (x2*y2 + x2 * y2));
	}

	double getValueB(double x1, double y1, double z1, double x2, double y2, double z2) {
		return 2 * (a*x1*x2 + b * y1*y2 + c * z1*z2 + f * (y1*z2 + y2 * z1) + g * (z1*x2 + z2 * x1) + h * (x1*y2 + x2 * y1) + p * (x2)+q * (y2)+r * (z2));
	}

	double getValueC(double x1, double y1, double z1, double x2, double y2, double z2) {
		return (a*x1*x1 + b * y1*y1 + c * z1*z1 + f * (y1*z1 + y1 * z1) + g * (z1*x1 + z1 * x1) + h * (x1*y1 + x1 * y1) +
			p * (x1 + x1) + q * (y1 + y1) + r * (z1 + z1) + d);
	}

	Vector getNormal(Point pt) {
		double i = a * (pt.x) + g * (pt.z) + h * (pt.y) + p; // factor of 2 is removed as normalization will be done
		double j = b * (pt.y) + f * (pt.z) + h * (pt.x) + q;
		double k = c * (pt.z) + f * (pt.y) + g * (pt.x) + r;
		double norm = sqrt(i * i + j * j + k * k);
		Vector normal(i / norm, j / norm, k / norm);
		return getUnitVector(multiplyMatrixVector(getNormalTransformationMatrix(), normal));
	}

	pair<int, int> getImageCoordinates(Point p) override {
		// texture is not handled for general quadric
		return { 0,0 };
	}

public:
	Quadric(double d1, double d2, double d3, double d4, double d5, double d6,
		double d7, double d8, double d9, double d10, Material m, vector<vector<double> > t, Point c)
		: Object(m), a(d1), b(d2), c(d3), f(d4), g(d5), h(d6), p(d7), q(d8), r(d9), d(d10), center(c)
	{
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