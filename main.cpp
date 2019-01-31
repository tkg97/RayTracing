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
			vector<Ray> rays = sceneObject.getRayFromViewer(i,j);
			vector<double> temp = { 0,0,0 };
			for (int k = 0; k < rays.size(); k++) {
				rgb = sceneObject.getIllumination(rays[k], -1, rdepth);
				for (int l = 0; l < 3; l++) {
					temp[l] += rgb[l];
				}
			}
			pixels[j*w + i].r = (temp[0]/5 <= 1) ? temp[0]/5 : 1;
			pixels[j*w + i].g = (temp[1]/5 <= 1) ? temp[1]/5 : 1;
			pixels[j*w + i].b = (temp[2]/5 <= 1) ? temp[2]/5 : 1;
        } 
    }
	
    savebmp("rayTrace.bmp", w, h, 72, pixels);

    return 0;
}
