#pragma once

class Material {
	vector<double> ambientCoefficient;
	vector<double> diffuseCoefficient;
	vector<double> specularCoefficient;
	double refractiveIndex;
	double phongExponent;
	double reflectionConstant;
	double refractionConstant;
	vector<double> textureImage; // bgr image
	int imageHeight, imageWidth;
public:
	Material(vector<double> a, vector<double> d, vector<double> s, double r, double p, double c1, double c2, string texturePath = "")
		: ambientCoefficient(a), diffuseCoefficient(d), specularCoefficient(s),
		refractiveIndex(r), phongExponent(p), reflectionConstant(c1), refractionConstant(c2)
	{
		if (texturePath != "") {
			unsigned char header[54]; // Each BMP file begins by a 54-bytes header
			unsigned int dataPos;     // Position in the file where the actual data begins
			unsigned int width, height;
			unsigned int imageSize;   // = width*height*3
			// Actual RGB data
			unsigned char * data;
			const char* tPath = texturePath.c_str();
			FILE * file = fopen(tPath, "rb");
			if (!file) { printf("Image could not be opened\n"); }
			if (fread(header, 1, 54, file) != 54) { // If not 54 bytes read : problem
				printf("Not a correct BMP file\n");
			}
			if (header[0] != 'B' || header[1] != 'M') {
				printf("Not a correct BMP file\n");
			}
			dataPos = *(int*)&(header[0x0A]);
			imageSize = *(int*)&(header[0x22]);
			width = *(int*)&(header[0x12]);
			height = *(int*)&(header[0x16]);
			if (imageSize == 0)    imageSize = width * height * 3; // 3 : one byte for each Red, Green and Blue component
			if (dataPos == 0)      dataPos = 54;
			data = new unsigned char[imageSize];

			// Read the actual data from the file into the buffer
			fread(data, 1, imageSize, file);

			//Everything is in memory now, the file can be closed
			fclose(file);
			imageHeight = height;
			imageWidth = width;
			for (int i = 0; i < imageSize; i++) {
				textureImage.push_back(((int)data[i]) / 255.0);
			}
		}
	}
	vector<double> getAmbientCoefficient() {
		return ambientCoefficient;
	}
	vector<double> getDiffusionCoefficeint() {
		return diffuseCoefficient;
	}
	vector<double> getSpecularCoefficient() {
		return specularCoefficient;
	}
	vector<double> getTextureCoefficeint(pair<int, int> imageMapUV) {
		// p will be provided in coordinate system of object
		if (!isTextureDefined()) return { 0,0,0 };
		else {
			if (imageMapUV.first >= imageWidth || imageMapUV.second > imageHeight || imageMapUV.second <= 0 || imageMapUV.first < 0) {
				return { 0,0,0 };
			}
			int g = 3 * (imageHeight - imageMapUV.second)*imageWidth + 3 * (imageMapUV.first) + 1;
			int r = 3 * (imageHeight - imageMapUV.second)*imageWidth + 3 * (imageMapUV.first) + 2;
			int b = 3 * (imageHeight - imageMapUV.second)*imageWidth + 3 * (imageMapUV.first) + 0;
			return { textureImage[r], textureImage[g], textureImage[b] };
		}
	}
	bool isTextureDefined() {
		return (!(textureImage.size() == 0));
	}
	double getRefractiveIndex() {
		return refractiveIndex;
	}
	double getPhongExponent() {
		return phongExponent;
	}
	double getReflectionConstant() {
		return reflectionConstant;
	}
	double getRefractionConstant() {
		return refractionConstant;
	}
	int getImageHeight() {
		return imageHeight;
	}
	int getImageWidth() {
		return imageWidth;
	}
};