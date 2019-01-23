#pragma once

#include <bits/stdc++.h>
#include "object.h"
#include "light.h"
#include "camera.h"
using namespace std;

class Scene{
        vector<Object> objects;
        vector<LightSource> lightSources;
        Viewer eye;
        double ambientLight;
    public:
        Scene(vector<Object>& o, vector<LightSource>& l, double a, Viewer& v) : objects(o), lightSources(l), ambientLight(a), eye(v){

        }
        vector<Object> getObjectList(){
            return objects;
        }
        vector<LightSource> getLightSourceList(){
            return lightSources;
        }
        double getAmbientLight(){
            return ambientLight;
        }
}