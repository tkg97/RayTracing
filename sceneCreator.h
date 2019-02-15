#include "json/json.h"
#include "object.h" 
#include "light.h"
#include "scene.h"
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

Point parseJsonArrayIntoPoint(Json::Value& pointLocation) noexcept(false){
	if (pointLocation.isArray()==false || pointLocation.size()!=3) throw exception();
	Point p(pointLocation[0].asDouble(), pointLocation[1].asDouble(), pointLocation[2].asDouble());
	return p;
}

vector<double> parseJsonArrayIntoDoubleVector(Json::Value& vectorData, int n) noexcept(false) {
	if (vectorData.isArray() == false || vectorData.size() != n) throw exception();
	vector<double> result;
    for(int i=0;i<n;i++){
        result.push_back(vectorData[i].asDouble());
    }
	return result;
}

Material parseJsonIntoMaterial(Json::Value &materialData) noexcept(false) {
	int t = materialData.size();
	if(materialData.isMember("ambientCoefficients") && materialData.isMember("diffusionCoefficients") && materialData.isMember("phongExponent")
		&& materialData.isMember("reflectionConstant") && materialData.isMember("refractionConstant") && materialData.isMember("refractiveIndex")
		&& materialData.isMember("specularCoefficients")) {

		vector<double> ambientCoefficients = parseJsonArrayIntoDoubleVector(materialData.get("ambientCoefficients", {}), 3);
		vector<double> diffusionCoefficients = parseJsonArrayIntoDoubleVector(materialData.get("diffusionCoefficients", {}), 3);
		vector<double> specularCoefficients = parseJsonArrayIntoDoubleVector(materialData.get("specularCoefficients", {}), 3);
		double refractiveIndex = materialData.get("refractiveIndex", 0.0).asDouble();
		double phongExponent = materialData.get("phongExponent", 1.0).asDouble();
		double reflectionConstant = materialData.get("reflectionConstant", 0.0).asDouble();
		double refractionConstant = materialData.get("refractionConstant", 0.0).asDouble();
		string texturePath = materialData.get("texture", "").asString();
		Material m(ambientCoefficients, diffusionCoefficients, specularCoefficients, refractiveIndex, phongExponent, reflectionConstant, refractionConstant, texturePath);
		return m;
	}
	else throw exception();
}

vector<Point> parseJsonIntoCoordinates(Json::Value &data, int n) noexcept(false){
	if (data.isArray() == false || data.size() != n) throw exception();
	vector<Point> coordinates;
	for (int i = 0;i < n;i++) {
		coordinates.push_back(parseJsonArrayIntoPoint(data[i]));
	}
	return coordinates;
}

vector<vector<double>> parseJsonArrayIntoMatrix(Json::Value &data) noexcept(false){
    if(data.isArray()==false || data.size()!=4) throw exception();
    vector< vector<double> > result;
    for(int i=0;i<4;i++){
        result.push_back(parseJsonArrayIntoDoubleVector(data[i], 4));
    }
    return result;
}

Scene getSceneObject(){
	vector<Object*> objectList;
	vector<LightSource> lightSourceList;
    ifstream f("inputFiles/Mirror/objectInput.json");
    Json::Reader reader;
    Json::Value root;
    if(reader.parse(f, root)){
      Json::Value::iterator itr = root.begin();
      while(itr!=root.end()){
		  try {
			  string typeOfObject = itr->get("type", "").asString();
			  if (typeOfObject == "sphere") {
				  Json::Value parameters = itr->get("parameters", {});
				  if (parameters.size()==0) throw exception();
				  double radius = parameters.get("radius", 0.0).asDouble();
				  if (radius <= 0) throw exception();
				  Json::Value centerLocation = parameters.get("center", {});
				  Point center = parseJsonArrayIntoPoint(centerLocation);
				  Json::Value materialDetails = itr->get("material", {});
				  Material objectMaterial = parseJsonIntoMaterial(materialDetails);
				  Json::Value affineTransformation = itr->get("affineTransformation", {});
				  vector<vector<double>> transformationMatrix = parseJsonArrayIntoMatrix(affineTransformation);
				  objectList.push_back(new Sphere(radius, center, objectMaterial, transformationMatrix));
			  }
			  else if(typeOfObject == "box"){
				  Json::Value parameters = itr->get("parameters", {});
				  if (parameters.size() == 0) throw exception();
				  vector<Point> coordinates = parseJsonIntoCoordinates(parameters.get("coordinates", {}), 8 /*number of vertices in box*/);
				  Point centerLocation = parseJsonArrayIntoPoint(parameters.get("center", {}));
				  Json::Value materialDetails = itr->get("material", {});
				  Material objectMaterial = parseJsonIntoMaterial(materialDetails);
				  Json::Value affineTransformation = itr->get("affineTransformation", {});
				  vector<vector<double>> transformationMatrix = parseJsonArrayIntoMatrix(affineTransformation);
				  objectList.push_back(new Box(coordinates, objectMaterial, transformationMatrix, centerLocation));
			  }
			  else if (typeOfObject == "polygon") {
				  Json::Value parameters = itr->get("parameters", {});
				  if (parameters.size() == 0) throw exception();
				  int vertexCount = parameters.get("vertexCount", 0).asInt();
				  if (vertexCount < 3) throw exception();// mininum triangle with three vertices
				  vector<Point> coordinates = parseJsonIntoCoordinates(parameters.get("coordinates", {}), vertexCount);
				  Point centerLocation = parseJsonArrayIntoPoint(parameters.get("center", {}));
				  Json::Value materialDetails = itr->get("material", {});
				  Material objectMaterial = parseJsonIntoMaterial(materialDetails);
				  Json::Value affineTransformation = itr->get("affineTransformation", {});
				  vector<vector<double>> transformationMatrix = parseJsonArrayIntoMatrix(affineTransformation);
				  objectList.push_back(new Polygon(vertexCount, coordinates, objectMaterial, transformationMatrix, centerLocation));
			  }
			  else if (typeOfObject == "quadric") {
				  Json::Value parameters = itr->get("parameters", {});
				  if (parameters.size() == 0) throw exception(); // actually zero quadric, not considering this case//
				  double a = parameters.get("a", 0).asDouble();
				  double b = parameters.get("b", 0).asDouble();
				  double c = parameters.get("c", 0).asDouble();
				  double f = parameters.get("f", 0).asDouble();
				  double g = parameters.get("g", 0).asDouble();
				  double h = parameters.get("h", 0).asDouble();
				  double p = parameters.get("p", 0).asDouble();
				  double q = parameters.get("q", 0).asDouble();
				  double r = parameters.get("r", 0).asDouble();
				  double d = parameters.get("d", 0).asDouble();
				  Point centerLocation = parseJsonArrayIntoPoint(parameters.get("center", {}));
				  Json::Value materialDetails = itr->get("material", {});
				  Material objectMaterial = parseJsonIntoMaterial(materialDetails);
				  Json::Value affineTransformation = itr->get("affineTransformation", {});
				  vector<vector<double>> transformationMatrix = parseJsonArrayIntoMatrix(affineTransformation);
				  objectList.push_back(new Quadric(a, b, c, f, g, h, p, q, r, d, objectMaterial, transformationMatrix, centerLocation));
			  }
			  else {
				  throw exception();
			  }
		  }
		  catch (exception e) {
			  cerr << "invalid input file" << endl;
			  exit(0);
		  }
		  itr++;
      }
    }
	else {
		cerr << "object input file couldn't be opened" << endl;
		exit(0);
	}
	f.close();

	f.open("inputFiles/Mirror/lightInput.json");
	if (reader.parse(f, root)) {
		Json::Value::iterator itr = root.begin();
		while (itr != root.end()) {
			Point location = parseJsonArrayIntoPoint(itr->get("location", {}));
			vector<double> intensity = parseJsonArrayIntoDoubleVector(itr->get("intensity", {}), 3);
			LightSource l(location, intensity);
			lightSourceList.push_back(l);
			itr++;
		}
	}
	else {
		cerr << "light sources input file couldn't be opened" << endl;
		exit(0);
	}
	f.close();

    f.open("inputFiles/Mirror/viewerInput.json");
    if(reader.parse(f,root)){
        int width = root.get("width", 0.0).asInt();
        int height = root.get("height", 0.0).asInt();
        if(width<=0 || height <=0) throw exception();
        double angle = root.get("angle", 0.25).asDouble(); // in terms of multiple of M_PI // default set to 45 degrees
        vector<vector<double>> transformationMatrix = parseJsonArrayIntoMatrix(root.get("transformation", {}));
        vector<double> ambientLight = parseJsonArrayIntoDoubleVector(root.get("ambience", {}),3);
        int recursionDepth = root.get("depth",0).asInt();
        Viewer v(angle, width, height, transformationMatrix);
        Scene s(objectList, lightSourceList,ambientLight,v,recursionDepth);
        return s;
    }
	else {
		cerr << "viewer input file couldn't be opened" << endl;
		exit(0);
	}
	f.close();
}