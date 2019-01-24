#pragma once

#include <vector>
#include "object.h"
#include "light.h"
#include "camera.h"
using namespace std;

class Scene{
        vector<Object> objects;
        vector<LightSource> lightSources;
        Viewer eye;
        vector<double> ambientLight;
    public:
        Scene(vector<Object>& o, vector<LightSource>& l, vector<double> a, Viewer& v) : objects(o), lightSources(l), ambientLight(a), eye(v){

        }
        vector<Object> getObjectList(){
            return objects;
        }
        vector<LightSource> getLightSourceList(){
            return lightSources;
        }
        vector<double> getAmbientLight(){
            return ambientLight;
        }
        vector<double> getIlumination(Ray& r, int objID){
            // objID is the index of the object from which the ray is emerged
            // first check whether this ray intersects with any of the object or not
            // if no intersection than, [0,0,0] is trivially returned;
            IntersectionPoint* minIntersectionPoint = nullptr;
            double minIntersectionParameter = DBL_MAX;
            int objectIndex = -1; // Will correspond to the first intersected object
            for(int i=0;i<objects.size();i++){
                IntersectionPoint* p;
                if(i==objID){
                    p = objects[i].getIntersection(r, 0.001); // t=0 intersection is not required
                }
                else p = objects[i].getIntersection(r, 0);
                if(p!=nullptr && p->getRayParameter() < minIntersectionParameter){
                    minIntersectionParameter = p->getRayParameter();
                    minIntersectionPoint = p;
                    objectIndex = i;
                }
            }
            if(minIntersectionPoint==nullptr){
                return vector<double>(3,0);
            }
            else{
                // now we have point of intersection, calculate the illumination
                // Let us first check for shadows
                vector<double> illumination(3,0);
                for(int i=0;i<lightSources.size();i++){
                    double shadowParameter = 0; // 0 means no shadow, 1 means complete shadow
                    Vector direction = getSubtractionVector(minIntersectionPoint->getLocation(), lightSources[i].getLocation());
                    Ray shadowRay(minIntersectionPoint->getLocation(), direction);
					for (int j = 0; j< objects.size();j++) {
						IntersectionPoint* p;
                        if(j==objectIndex){
                            p = objects[j].getIntersection(shadowRay, 0.001);
                        }
                        else p = objects[j].getIntersection(shadowRay, 0);
						//TODO : currently shadow parameter is 0/1, update it with heuristic of opaqueness
						if (p != nullptr) {
							shadowParameter = 1;
						}
					}
                    if(shadowParameter!=1){
                        Vector unitLightDirection = getUnitVector(direction);
                        double cosineLightNormal = dotProduct(unitLightDirection,minIntersectionPoint->getNormal());
                        vector<double> diffusionIllumination = multiplyVectorDouble(cosineLightNormal,
                            multiplyVectorsPointwise(objects[objectIndex].getDiffusionCoefficeint(), lightSources[i].getIntensity()));
                        Vector unitViewDirection = getUnitVector(getSubtractionVector(minIntersectionPoint->getLocation(), eye.getEyeLocation()));
                        double cosineLightView = dotProduct(unitLightDirection, unitViewDirection);
                        double cosineNormalView = dotProduct(minIntersectionPoint->getNormal(), unitViewDirection);
                        double cosineReflectionView = 2*(cosineLightNormal)*(cosineNormalView) - (cosineLightView); 
                        vector<double> specularIllumination = multiplyVectorDouble(pow(cosineReflectionView, objects[objectIndex].getPhongExponent()),
                            multiplyVectorsPointwise(objects[objectIndex].getSpecularCoefficient(), lightSources[i].getIntensity()));
                        // update total illumination
                        for(int i=0;i<3;i++){
                            illumination[i] += diffusionIllumination[i] + specularIllumination[i];
                        }
                    } 
                }
				// now only ambient illumination is to be handled
                vector<double> ambientIllumination = multiplyVectorsPointwise(objects[objectIndex].getAmbientCoefficient(), ambientLight);
                // update total illumination
                for(int i=0;i<3;i++){
                    illumination[i] += ambientIllumination[i];
                }
                // now phong illumination is done at this point
                // reflection and refraction has to be handled now
                
            }
        }
};