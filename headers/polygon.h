#pragma once

class Polygon : public Object {
	int vertexCount; //vertexCount
	vector< Point > coordinates;
	Vector normal;
	Point center;


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
		return getUnitVector(multiplyMatrixVector(getNormalTransformationMatrix(), normal));
	}
	pair<int, int> getImageCoordinates(Point p) override {
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
		int u, v;
		if (param1[0] <= 1.0001 && param1[0] >= -0.0001 && param1[0] <= 1.0001 && param1[1] >= -0.0001 && param1[2] <= 1.0001 && param1[2] >= -0.0001 && (param1[0] + param1[1] + param1[2]) <= 1.0001)
		{
			u = (int)((param1[1] + param1[2])*objectMat.getImageWidth());
			v = (int)(param1[2] * objectMat.getImageHeight());
		}
		else {
			u = (int)((param2[0])*objectMat.getImageWidth());
			v = (int)((param2[0] + param2[1]) * objectMat.getImageHeight());
		}
		return { u,v };

	}

public:
	Polygon(int n, vector<Point> v, Material m, vector<vector<double>> t, Point c) :
		Object(m, true), vertexCount(n), coordinates(v), normal(0, 0, 0), center(c) {
		Point p1 = coordinates[0];
		Point p2 = coordinates[1];
		Point p3 = coordinates[2];

		normal = getUnitVector(crossProduct(p1, p2, p3));
		vector<vector<double>> translationMat = formTranslationMatrix(center.x, center.y, center.z);
		setTransformations(t, translationMat);
	}

	vector<Point> getCoordinates() {
		return coordinates;
	}
	// get the intersection of ray with the polygon
	IntersectionPoint* getIntersection(Ray r1, double minThreshold) {
		r1 = getTransformedRay(r1, getInverseVertexTransformationMatrix(), getRayTransformation());
		Point r0 = r1.getSource();
		Point p0 = coordinates[0];
		Vector rd = r1.getDirection();
		double a1 = r0.x - p0.x;
		double a2 = r0.y - p0.y;
		double a3 = r0.z - p0.z;
		Vector v(a1, a2, a3);
		if (abs(dotProduct(normal, rd)) <= 0.00001) {
			return nullptr;
		}
		double t = -(dotProduct(normal, v)) / (dotProduct(normal, rd));
		if (t < minThreshold) return nullptr;

		Point intersection = (addPointVector(r0, multiplyVectorDouble(t, rd)));
		if (isContained(intersection)) {
			return (new IntersectionPoint(multiplyMatrixVector(getVertexTransformation(), intersection), getNormal()));
		}
		else
			return nullptr;
	}
};