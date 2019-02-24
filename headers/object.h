#pragma once
#define _USE_MATH_DEFINES
#include <vector>
#include <iostream>
#include <algorithm>
#include <float.h>
#include "point.h"
#include "matrixAlgebra.h"
#include "material.h"
using namespace std;

class Object {
	// an abstract parent class for all the objects
	Material objectMaterial;
	bool planarity;
	vector<vector<double>> vertexTransformation;
	vector<vector<double>> inverseVertexTransformation;
	vector<vector<double>> normalTransformation;
	vector<vector<double>> rayTransformation;
	virtual pair<int, int> getImageCoordinates(Point p) = 0;
public:
	Object(Material m, bool planar = false) : objectMaterial(m), planarity(planar) {}
	vector<double> getAmbientCoefficient() {
		return objectMaterial.getAmbientCoefficient();
	}
	vector<double> getDiffusionCoefficeint(Point p) {
		if (objectMaterial.isTextureDefined()) {
			p = multiplyMatrixVector(getInverseVertexTransformationMatrix(), p);
			pair<int, int> imageMapUV = getImageCoordinates(p);
			return objectMaterial.getTextureCoefficeint(imageMapUV);
		}
		return objectMaterial.getDiffusionCoefficeint();
	}
	vector<double> getSpecularCoefficient() {
		return objectMaterial.getSpecularCoefficient();
	}
	double getRefractiveIndex() {
		return objectMaterial.getRefractiveIndex();
	}
	double getPhongExponent() {
		return objectMaterial.getPhongExponent();
	}
	double getReflectionConstant() {
		return objectMaterial.getReflectionConstant();
	}
	double getRefractionConstant() {
		return objectMaterial.getRefractionConstant();
	}
	bool isPlanar() {
		return planarity;
	}
	vector<vector<double>> getVertexTransformation() {
		return vertexTransformation;
	}
	vector<vector<double>> getInverseVertexTransformationMatrix() {
		return inverseVertexTransformation;
	}
	vector<vector<double>> getNormalTransformationMatrix() {
		return normalTransformation;
	}
	vector<vector<double>> getRayTransformation() {
		return rayTransformation;
	}
	Material getObjectMaterial() {
		return objectMaterial;
	}
	virtual IntersectionPoint* getIntersection(Ray r, double minThreshold) = 0;
	// minThreshold will help me handle the case for t==0 avoidance
protected:
	void setTransformations(const vector<vector<double>> &t, vector<vector<double>> &trans) {
		rayTransformation = getMatrixInverse(t);
		vertexTransformation = multiplyMatrices(getMatrixInverse(trans), multiplyMatrices(t, trans));
		inverseVertexTransformation = multiplyMatrices(getMatrixInverse(trans), multiplyMatrices(rayTransformation, trans));
		normalTransformation = getMatrixTranspose(rayTransformation);
	}
};