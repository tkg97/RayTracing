#pragma once
#include <vector>
#include <cfloat>
#include "object.h"
#include "light.h"
#include <cmath>
#include <limits.h>
#include <limits>
#include "camera.h"
using namespace std;

class Scene{
        vector<Object*> objects;
        vector<LightSource> lightSources;
        Viewer eye;
        vector<double> ambientLight;
        int recursionDepth;
		double scale;
		double imageAspectRatio;

        Vector getTransmissionVector(Vector i, Vector n, double rIndex1, double rIndex2){
            double c = abs(dotProduct(i,n));
            double r = rIndex1/rIndex2;
            Vector direction = multiplyVectorDouble(r, i) + multiplyVectorDouble(-sqrt(1 - r*r*(1 - c*c)) + r*c, n);
            return getUnitVector(direction);
        }
    public:
        Scene(vector<Object*> o, vector<LightSource> l, vector<double> a, Viewer v, int r) : objects(o), lightSources(l), ambientLight(a), eye(v), recursionDepth(r){
			scale = tan(eye.getAngle()*M_PI);
			imageAspectRatio = (eye.getWidth()*1.0) / (eye.getHeight());
		}
        vector<Object*> getObjectList(){
            return objects;
        }
        vector<LightSource> getLightSourceList(){
            return lightSources;
        }
        vector<double> getAmbientLight(){
            return ambientLight;
        }
        vector<double> getIllumination(Ray& r, int objID, int depth){
            // objID is the index of the object from which the ray is emerged
            // if depth becomes zero, reflection and transmission is not considered
            // toEnter differentiates whether it is an entering or exiting ray
            // first check whether this ray intersects with any of the object or not
            // if no intersection than, [0,0,0] is trivially returned;
            
            IntersectionPoint* minIntersectionPoint = nullptr;
            double minIntersectionParameter = DBL_MAX;
            int objectIndex = -1; // Will correspond to the first intersected object
            for(int i=0;i<objects.size();i++){
                IntersectionPoint* p;
                if(i==objID){
                    p = objects[i]->getIntersection(r, 0.001); // t=0 intersection is not required
                }
                else p = objects[i]->getIntersection(r, 0);
                if(p!=nullptr && p->getRayParameter() < minIntersectionParameter){
                    minIntersectionParameter = p->getRayParameter();
                    minIntersectionPoint = p;
                    objectIndex = i;
                }
            }
            if(minIntersectionPoint==nullptr){
                vector<double> illumination(3,0);
				return illumination;
            }
            else{
                // now we have point of intersection, calculate the illumination
                // Let us first check for shadows
				if (dotProduct(r.getDirection(), minIntersectionPoint->getNormal()) >= 0 && typeid(*(objects[objectIndex])) != typeid(Polygon)) {
					// means the ray is exiting
					Vector transmissionDirection = getTransmissionVector(r.getDirection(), minIntersectionPoint->getNormal(), objects[objectIndex]->getRefractiveIndex(), 1);
					Ray transmissionRay(minIntersectionPoint->getLocation(), transmissionDirection);
					vector<double> illumination = getIllumination(transmissionRay, objectIndex, depth);
					return illumination;
				}

                vector<double> illumination(3,0);
                for(int i=0;i<lightSources.size();i++){
                    int shadowParameter = 0; // 0 means no shadow, 1 means complete shadow
                    Vector direction = getSubtractionVector(minIntersectionPoint->getLocation(), lightSources[i].getLocation());
                    double lightSourceDistance = sqrt((direction.i * direction.i) + (direction.j * direction.j) + (direction.k * direction.k));
					Ray shadowRay(minIntersectionPoint->getLocation(), direction);
					double cosineRayNormal = dotProduct(r.getDirection(), minIntersectionPoint->getNormal());
					double cosineLightNormal = dotProduct(shadowRay.getDirection(), minIntersectionPoint->getNormal());
					if ((cosineLightNormal >= 0 && cosineRayNormal >= 0) || (cosineLightNormal <= 0 && cosineRayNormal <= 0)) shadowParameter = 1;
					for (int j = 0; j< objects.size();j++) {
						IntersectionPoint* p;
                        if(j==objectIndex){
                            p = objects[j]->getIntersection(shadowRay, 0.001);
                        }
                        else p = objects[j]->getIntersection(shadowRay, 0);
						//TODO : currently shadow parameter is 0/1, update it with heuristic of opaqueness
						if (p != nullptr && p->getRayParameter()<=lightSourceDistance) {
							shadowParameter = 1;
						}
					}
                    if(shadowParameter!=1){
                        Vector unitLightDirection = getUnitVector(direction); //shadowRay.getDirection()
                        double cosineLightNormal = abs(dotProduct(unitLightDirection,minIntersectionPoint->getNormal()));
                        vector<double> diffusionIllumination = multiplyVectorDouble(cosineLightNormal,
                            multiplyVectorsPointwise(objects[objectIndex]->getDiffusionCoefficeint(), lightSources[i].getIntensity()));
						Vector unitViewDirection = multiplyVectorDouble(-1, r.getDirection());
						double cosineLightView = dotProduct(unitLightDirection, unitViewDirection);
                        double cosineNormalView = abs(dotProduct(minIntersectionPoint->getNormal(), unitViewDirection));
                        double cosineReflectionView = 2*(cosineLightNormal)*(cosineNormalView) - (cosineLightView); 
                        vector<double> specularIllumination = multiplyVectorDouble(pow(cosineReflectionView, objects[objectIndex]->getPhongExponent()),
                            multiplyVectorsPointwise(objects[objectIndex]->getSpecularCoefficient(), lightSources[i].getIntensity()));
                        // update total illumination
                        for(int i=0;i<3;i++){
                            illumination[i] += (diffusionIllumination[i] + specularIllumination[i]);
                        }
                    } 
                }
				// now only ambient illumination is to be handled
                if(depth==recursionDepth){
                    vector<double> ambientIllumination = multiplyVectorsPointwise(objects[objectIndex]->getAmbientCoefficient(), ambientLight);
                    // update total illumination
                    for(int i=0;i<3;i++){
                        illumination[i] += ambientIllumination[i];
                    }
                }
                // now phong illumination is done at this point
                // reflection and refraction has to be handled now
                
                if(depth !=0){
                    //Reflection
                    if(objects[objectIndex]->getReflectionConstant()>=0.0001){
                        double cosineIncidentNormal = dotProduct(r.getDirection(), minIntersectionPoint->getNormal());
                        Vector reflectionDirection(0,0,0);
                        reflectionDirection = r.getDirection() + multiplyVectorDouble(-2*cosineIncidentNormal, minIntersectionPoint->getNormal());
                        Ray reflectionRay(minIntersectionPoint->getLocation(), reflectionDirection);

                        vector<double> reflectionIllumination = getIllumination(reflectionRay, objectIndex, depth-1);
                        // update total illumination
                        for(int i=0;i<3;i++){
                            illumination[i] += objects[objectIndex]->getReflectionConstant()*reflectionIllumination[i];
                        }
                    }

                    //Refraction
                    if(objects[objectIndex]->getRefractionConstant()>=0.0001){
                        if(objects[objectIndex]->isPlanar()){
                            // if object is planar, than no refraction is considered (case of polygon), simple transmission
                            Ray transmissionRay(minIntersectionPoint->getLocation(), r.getDirection());
                            vector<double> transmissionIllumination = getIllumination(transmissionRay, objectIndex, depth-1);
                            // update total illumination
                            for(int i=0;i<3;i++){
                                illumination[i] += objects[objectIndex]->getRefractionConstant()*transmissionIllumination[i];
                            }
                        }
                        else{
                            Vector transmissionDirection = getTransmissionVector(r.getDirection(), minIntersectionPoint->getNormal(), 1, objects[objectIndex]->getRefractiveIndex());
                            Ray transmissionRay(minIntersectionPoint->getLocation(), transmissionDirection);
                            vector<double> transmissionIllumination = getIllumination(transmissionRay, objectIndex, depth-1);
                            // update total illumination
                            for(int i=0;i<3;i++){
                                illumination[i] += objects[objectIndex]->getRefractionConstant()*transmissionIllumination[i];
                            }
                        }
                    }
                }
                return illumination;
            }
        }
        
		double getViewingAngle(){
            return eye.getAngle();
        }
        vector<vector<double> > getTransformationMatrix(){
            return eye.getTransformationMatrix();
        }
        int getRecursionDepth(){
            return recursionDepth;
        }
        int getWindowWidth(){
            return eye.getWidth();
        }
        int getWindowHeight(){
            return eye.getHeight();
        }
        vector<Ray> getRayFromViewer(int i, int j){
            double x = (2 * ((i + 0.5) / (eye.getWidth()*1.0)) - 1) * imageAspectRatio * scale; 
            double y = (2 * ((j + 0.5) / (eye.getHeight()*1.0)) - 1) * scale;
			double x1 = x + scale / (eye.getHeight()*1.0);
			double x2 = x - scale / (eye.getHeight()*1.0);
			double y1 = y + scale / (eye.getHeight()*1.0);
			double y2 = y - scale / (eye.getHeight()*1.0);
			Vector dir1 = eye.getRayDirection(x, y);
			Vector dir2 = eye.getRayDirection(x1, y1);
			Vector dir3 = eye.getRayDirection(x1, y2);
			Vector dir4 = eye.getRayDirection(x2, y1);
			Vector dir5 = eye.getRayDirection(x2, y2);
            Ray r1(eye.getEyeLocation() , dir1);
			Ray r2(eye.getEyeLocation(), dir2);
			Ray r3(eye.getEyeLocation(), dir3);
			Ray r4(eye.getEyeLocation(), dir4);
			Ray r5(eye.getEyeLocation(), dir5);
			vector<Ray> rays;
			rays.push_back(r1);
			rays.push_back(r2);
			rays.push_back(r3);
			rays.push_back(r4);
			rays.push_back(r5);
            return rays;
        }
};