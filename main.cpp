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

int render(const std::vector<float>& originalRayData, const std::vector<float>& shadowRayData, 
	const std::vector<float>& reflectedRayData, const std::vector<float>& refractedRayData);

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
				if (k == 0) {
					rgb = sceneObject.getIllumination(rays[k], -1, rdepth, j, i, rayType::Original);
				}
				else rgb = sceneObject.getIllumination(rays[k], -1, rdepth, -1, -1, rayType::Original); // so that only one ray is considered for opengl visualisation
				for (int l = 0; l < 3; l++) {
					temp[l] += rgb[l];
				}
			}
			pixels[j*w + i].r = (temp[0] / rays.size() <= 1) ? temp[0] / rays.size() : 1;
			pixels[j*w + i].g = (temp[1] / rays.size() <= 1) ? temp[1] / rays.size() : 1;
			pixels[j*w + i].b = (temp[2] / rays.size() <= 1) ? temp[2] / rays.size() : 1;
        } 
    }
	
    savebmp("inputFiles/Opengl/rayTrace.bmp", w, h, 72, pixels);

	cout << "Please enter the pixel for which you want to see the simulation" << endl;
	cout << "Enter the two integers i and j" << endl;
	cout << "Enter -1 to stop entering the values" << endl;
	while (1) {
		int pixelI, pixelJ;
		cin >> pixelI;
		if (pixelI == -1) break;
		cin >> pixelJ;
		if (pixelJ == -1) break;

		vector<float> originalPixelData, shadowPixelData, reflectedPixelData, refractedPixelData;
		sceneObject.getRequiredPixelData(pixelI, pixelJ, originalPixelData, shadowPixelData, reflectedPixelData, refractedPixelData);
		render(originalPixelData, shadowPixelData, reflectedPixelData, refractedPixelData);
	}
    return 0;
}
