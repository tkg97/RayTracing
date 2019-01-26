#include<vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "scene.h"
#include "point.h"
#include "camera.h"
#include "object.h"
#include "light.h"
using namespace std;

struct RGBType {
	double r;
	double g;
	double b;
};


//producing a bitmap image
void savebmp(const char *filename, int w, int h, int dpi, RGBType *data) {
	FILE *f;
	int k = w*h;
	int s = 4 * k;
	int filesize = 54 + s;

	double factor = 39.375;
	int m = static_cast<int>(factor);

	int ppm = dpi*m;

	unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0,0,0, 54,0,0,0 };
	unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0 };

	bmpfileheader[2] = (unsigned char)(filesize);
	bmpfileheader[3] = (unsigned char)(filesize >> 8);
	bmpfileheader[4] = (unsigned char)(filesize >> 16);
	bmpfileheader[5] = (unsigned char)(filesize >> 24);

	bmpinfoheader[4] = (unsigned char)(w);
	bmpinfoheader[5] = (unsigned char)(w >> 8);
	bmpinfoheader[6] = (unsigned char)(w >> 16);
	bmpinfoheader[7] = (unsigned char)(w >> 24);

	bmpinfoheader[8] = (unsigned char)(h);
	bmpinfoheader[9] = (unsigned char)(h >> 8);
	bmpinfoheader[10] = (unsigned char)(h >> 16);
	bmpinfoheader[11] = (unsigned char)(h >> 24);

	bmpinfoheader[21] = (unsigned char)(s);
	bmpinfoheader[22] = (unsigned char)(s >> 8);
	bmpinfoheader[23] = (unsigned char)(s >> 16);
	bmpinfoheader[24] = (unsigned char)(s >> 24);

	bmpinfoheader[25] = (unsigned char)(ppm);
	bmpinfoheader[26] = (unsigned char)(ppm >> 8);
	bmpinfoheader[27] = (unsigned char)(ppm >> 16);
	bmpinfoheader[28] = (unsigned char)(ppm >> 24);

	bmpinfoheader[29] = (unsigned char)(ppm);
	bmpinfoheader[30] = (unsigned char)(ppm >> 8);
	bmpinfoheader[31] = (unsigned char)(ppm >> 16);
	bmpinfoheader[32] = (unsigned char)(ppm >> 24);

	f = fopen(filename, "wb");

	fwrite(bmpfileheader, 1, 14, f);
	fwrite(bmpinfoheader, 1, 40, f);

	for (int i = 0; i < k; i++) {
		RGBType rgb = data[i];

		double red = data[i].r * 255;
		double green = (data[i].g) * 255;
		double blue = (data[i].b) * 255;

		unsigned char color[3] = { (unsigned char)floor(blue),(unsigned char)floor(green),(unsigned char)floor(red) };

		fwrite(color, 1, 3, f);
	}

	fclose(f);
}



// multiply a matrix(4*4) with an augmented point to get the point in the new coordinate system
Point multiplyMatrix(vector<vector<double>> cameraToWorld,Point p){
    double u1 = p.x*cameraToWorld[0][0] + p.y*cameraToWorld[1][0] + p.z*cameraToWorld[2][0] + cameraToWorld[3][0];
    double u2 = p.x*cameraToWorld[0][1] + p.y*cameraToWorld[1][1] + p.z*cameraToWorld[2][1] + cameraToWorld[3][1];
    double u3 = p.x*cameraToWorld[0][2] + p.y*cameraToWorld[1][2] + p.z*cameraToWorld[2][2] + cameraToWorld[3][2];
    Point r(u1,u2,u3);
    return r;
}
int main(){
    //creating scene
    vector<vector<double>> cameraToWorld = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,0}}; 
    Point o(0,0,0);
    Point origin = multiplyMatrix(cameraToWorld, o);
    int w = 1024; 
	int h = 768;
    double angle = 45; // in degrees
    Viewer v(origin,angle, w,h);
    vector<Object*> spheres;
    Point centre1(0.0, 0, -20);
    double rad1 = 4;
    Sphere s1(rad1 ,centre1, {0.2, 0.2, 0.2},{0.8, 0.8, 0.8},{0.1,0.1,0.1}, 5, 20 , 0.5, 0.1);
    Point centre2(0.0, 10, -20);
    double rad2 = 4;
    Sphere s2(rad2 ,centre2, {0.2, 0.2, 0.2},{0.8, 0.8, 0.8},{0.1,0.1,0.1}, 5, 20 , 0.5, 0.1);
    spheres.push_back(&s1);
    spheres.push_back(&s2);
    vector<LightSource> lightSources;
    Point Location1(0,20,-7);
    LightSource l1(Location1,{1,1,1});
	lightSources.push_back(l1);
    Scene s(spheres, lightSources, {1,1,1}, v );
    

    //creating the window
    RGBType *pixels = new RGBType[w*h];


    // rendering the scene
    double scale = tan((angle * M_PI)/180); 
    double imageAspectRatio = (w*1.0)/h; 
    for (int i = 0; i < w; ++i) { 
        for (int j = 0; j < h; ++j) { 
            float x = (2 * ((i + 0.5) / (w*1.0)) - 1) * imageAspectRatio * scale; 
            float y = (2 * ((j + 0.5) / (h*1.0)) - 1) * scale; 
            Point a(x,y,-1);
            Point b = multiplyMatrix(cameraToWorld ,a );
            Vector dir = getSubtractionVector(origin,b); 
            Ray r(origin , dir);
            vector<double> rgb= s.getIllumination(r, -1, 0, true);
            pixels[j*w + i].r = rgb[0];
			pixels[j*w + i].g = rgb[1];
			pixels[j*w + i].b = rgb[2];
        } 
    }
	// for (int i = 0; i < w; ++i) { 
    //     for (int j = 0; j < h; ++j) {
	// 		cout << "[" << pixels[j*w + i].r << " " << pixels[j*w + i].g << " " << pixels[j*w + i].b << "]" << " ";
	// 	}
	// 	cout << endl;
	// }
	cout << noIntersect << endl;
	cout << shadow << endl;
    savebmp("rayTrace.bmp", w, h, 72, pixels);

    return 0;
}


