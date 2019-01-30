#define _USE_MATH_DEFINES
#include<vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "scene.h"
#include "point.h"
#include "camera.h"
#include "object.h"
#include "light.h"
#include "sceneCreator.h"
#include "bmpImageHandler.h"
using namespace std;

int main(){
    //creating scene
    Scene sceneObject = getSceneObject(); // read input files to get scene object
    //creating the window
	int w = sceneObject.getWindowWidth();
	int h = sceneObject.getWindowHeight();
    RGBType *pixels = new RGBType[w*h];

    // rendering the scene

	int rdepth = sceneObject.getRecursionDepth();
	vector<double> rgb;

    for (int i = 0; i < w; ++i) {
        for (int j = 0; j < h; ++j) {
			Ray r = sceneObject.getRayFromViewer(i,j);
            rgb = sceneObject.getIllumination(r, -1, rdepth, true);
            pixels[j*w + i].r = (rgb[0]<=1) ? rgb[0] : 1;
            pixels[j*w + i].g = (rgb[1]<=1) ? rgb[1] : 1;
            pixels[j*w + i].b = (rgb[2]<=1) ? rgb[2] : 1;
        } 
    }
	cout << noIntersect << endl;
	cout << shadow << endl;
    savebmp("rayTrace.bmp", w, h, 72, pixels);

    return 0;
}

//TODO: light source into the object
