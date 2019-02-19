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

int render(const std::vector<float>& rayData);

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
				rgb = sceneObject.getIllumination(rays[k], -1, rdepth, j, i);
				for (int l = 0; l < 3; l++) {
					temp[l] += rgb[l];
				}
			}
			pixels[j*w + i].r = (temp[0] / rays.size() <= 1) ? temp[0] / rays.size() : 1;
			pixels[j*w + i].g = (temp[1] / rays.size() <= 1) ? temp[1] / rays.size() : 1;
			pixels[j*w + i].b = (temp[2] / rays.size() <= 1) ? temp[2] / rays.size() : 1;
        } 
    }
	
    savebmp("rayTrace.bmp", w, h, 72, pixels);

	cout << "Please enter the pixel for which you want to see the simulation" << endl;
	cout << "Enter the two integers i and j" << endl;
	cout << "Enter -1 to stop entering the values" << endl;
	while (1) {
		int pixelI, pixelJ;
		cin >> pixelI;
		if (pixelI == -1) break;
		cin >> pixelJ;
		if (pixelJ == -1) break;

		vector<float> reqDataPixel;
		reqDataPixel = sceneObject.getRequiredPixelData(pixelI, pixelJ);
		render(reqDataPixel);
	}
    return 0;
}
