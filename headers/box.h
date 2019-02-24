#pragma once
#include "polygon.h"

class Box : public Object {
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
	Point center;

	void storeAllPolygons(const Material &m, const vector<vector<double>> &t) {
		Polygon p1(4, { coordinates[0], coordinates[1], coordinates[2], coordinates[3] }, m, t, center);
		Polygon p2(4, { coordinates[7], coordinates[4], coordinates[0], coordinates[3] }, m, t, center);
		Polygon p3(4, { coordinates[6], coordinates[5], coordinates[4], coordinates[7] }, m, t, center);
		Polygon p4(4, { coordinates[2], coordinates[1], coordinates[5], coordinates[6] }, m, t, center);
		Polygon p5(4, { coordinates[5], coordinates[1], coordinates[0], coordinates[4] }, m, t, center);
		Polygon p6(4, { coordinates[7], coordinates[3], coordinates[2], coordinates[6] }, m, t, center);
		polygonFaces = vector<Polygon>({ p1, p2, p3, p4, p5, p6 });
	}

	pair<int, int> getImageCoordinates(Point p) override {
		int u, v;
		for (int i = 0; i < polygonFaces.size(); i++) {
			vector<Point> coordinates = polygonFaces[i].getCoordinates();
			double area1 = TriangleArea(p, coordinates[1], coordinates[2]);
			double area2 = TriangleArea(p, coordinates[0], coordinates[2]);
			double area3 = TriangleArea(p, coordinates[1], coordinates[0]);
			double Area = TriangleArea(coordinates[0], coordinates[1], coordinates[2]);

			double area4 = TriangleArea(p, coordinates[0], coordinates[3]);
			double area5 = TriangleArea(p, coordinates[0], coordinates[2]);
			double area6 = TriangleArea(p, coordinates[2], coordinates[3]);
			double Area1 = TriangleArea(coordinates[0], coordinates[2], coordinates[3]);

			vector<double> param1 = { area1 / Area, area2 / Area , area3 / Area };
			vector<double> param2 = { area4 / Area1, area5 / Area1 , area6 / Area1 };
			Material objectMat = getObjectMaterial();
			if (param1[0] <= 1.0001 && param1[0] >= -0.0001 && param1[0] <= 1.0001 && param1[1] >= -0.0001 && param1[2] <= 1.0001 && param1[2] >= -0.0001 && (param1[0] + param1[1] + param1[2]) <= 1.0001)
			{
				u = (int)((param1[1] + param1[2])*objectMat.getImageWidth());
				v = (int)(param1[2] * objectMat.getImageHeight());
			}
			else if (param2[0] <= 1.0001 && param2[0] >= -0.0001 && param2[0] <= 1.0001 && param2[1] >= -0.0001 && param2[2] <= 1.0001 && param2[2] >= -0.0001 && (param2[0] + param2[1] + param2[2]) <= 1.0001) {
				u = (int)((param2[0])*objectMat.getImageWidth());
				v = (int)((param2[0] + param2[1]) * objectMat.getImageHeight());
			}
		}
		return { u,v };
	}

public:
	Box(vector<Point> v, Material m, vector<vector<double> > t, Point c) : Object(m), coordinates(v), center(c) {
		vector<vector<double>> translationMat = formTranslationMatrix(center.x, center.y, center.z);
		setTransformations(t, translationMat);
		storeAllPolygons(m, t);
	}

	IntersectionPoint* getIntersection(Ray r, double minThreshold) {
		IntersectionPoint* minIntersectionPoint = nullptr;
		double minIntersectionParameter = DBL_MAX;
		for (int i = 0;i < 6;i++) {
			IntersectionPoint* intersection = polygonFaces[i].getIntersection(r, minThreshold);
			if (intersection != nullptr) {
				if (minIntersectionPoint == nullptr) {
					minIntersectionPoint = intersection;
					minIntersectionParameter = getNorm(getSubtractionVector(r.getSource(), intersection->getLocation()));
				}
				else if (minIntersectionParameter > getNorm(getSubtractionVector(r.getSource(), intersection->getLocation()))) {
					minIntersectionPoint = intersection;
					minIntersectionParameter = getNorm(getSubtractionVector(r.getSource(), intersection->getLocation()));

				}
			}
		}
		return minIntersectionPoint;
	}
};