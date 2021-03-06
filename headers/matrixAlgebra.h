#pragma once
#include "point.h"

double TriangleArea(Point x1, Point x2, Point x3)
{
	// Initialize area 
	double area = 0.0;

	Vector v1 = getSubtractionVector(x2, x1);
	Vector v2 = getSubtractionVector(x2, x3);
	// Calculate value of using dot product
	double ab = getNorm(v1);
	double ac = getNorm(v2);
	double cosTheta = dotProduct(v1, v2) / (ab * ac);
	double sinTheta = sqrt(1 - cosTheta * cosTheta);
	area = ab * ac* sinTheta;
	// Return absolute value 
	return abs(area / 2.0);
}

Point multiplyMatrixVector(const vector<vector<double>> &matrix, Point p) {
	// multiply a matrix(4*4) with an augmented point (input in 3d form only) to get the point in the new coordinate system
	double u1 = p.x*matrix[0][0] + p.y*matrix[1][0] + p.z*matrix[2][0] + 1 * matrix[3][0];
	double u2 = p.x*matrix[0][1] + p.y*matrix[1][1] + p.z*matrix[2][1] + 1 * matrix[3][1];
	double u3 = p.x*matrix[0][2] + p.y*matrix[1][2] + p.z*matrix[2][2] + 1 * matrix[3][2];
	Point r(u1, u2, u3);
	return r;
}

vector<double> multiplyMatrixPoint(const vector<vector<double>> &matrix, Point p) {
	// multiply a matrix(3*3) and a point to get the the paramters of the linear combination
	double u1 = p.x*matrix[0][0] + p.y*matrix[1][0] + p.z*matrix[2][0];
	double u2 = p.x*matrix[0][1] + p.y*matrix[1][1] + p.z*matrix[2][1];
	double u3 = p.x*matrix[0][2] + p.y*matrix[1][2] + p.z*matrix[2][2];
	vector<double> r = { u1, u2, u3 };
	return r;
}

Vector multiplyMatrixVector(const vector<vector<double>> &matrix, Vector v) {
	// multiply a matrix(4*4) with an augmented direction (input in 3d form only) to get the point in the new coordinate system
	double u1 = v.i*matrix[0][0] + v.j*matrix[1][0] + v.k*matrix[2][0] + 0 * matrix[3][0];
	double u2 = v.i*matrix[0][1] + v.j*matrix[1][1] + v.k*matrix[2][1] + 0 * matrix[3][1];
	double u3 = v.i*matrix[0][2] + v.j*matrix[1][2] + v.k*matrix[2][2] + 0 * matrix[3][2];
	Vector d(u1, u2, u3);
	return d;
}

vector<vector<double>> multiplyMatrices(const vector<vector<double>>& mat1, const vector<vector<double>> &mat2) {
	// both will be 4*4 in our code, but it is written for generic matrix multiplication;
	// both matrices are assumed to be compatible with each other
	int len1 = mat1.size();
	int len2 = mat1[0].size();
	int len3 = mat2[0].size();
	vector<vector<double>> resultMat;
	for (int i = 0;i < len1;i++) {
		vector<double> temprow;
		for (int k = 0;k < len3;k++) {
			double value = 0;
			for (int j = 0;j < len2;j++) {
				value += mat1[i][j] * mat2[j][k];
			}
			temprow.push_back(value);
		}
		resultMat.push_back(temprow);
	}
	return resultMat;
}

void getCofactor(vector<vector<double>> A, vector<vector<double>> &temp, int p, int q, int n)
{
	int i = 0, j = 0;
	double** temp1 = new double*[n];
	for (int k = 0; k < n; k++) {
		temp1[k] = new double[n];
	}
	// Looping for each element of the matrix 
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			//  Copying into temporary matrix only those element 
			//  which are not in given row and column 
			if (row != p && col != q)
			{
				temp1[i][j++] = A[row][col];

				// Row is filled, so increase row index and 
				// reset col index 
				if (j == n - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
	for (int i = 0; i < n; i++) {
		vector<double> temp2;
		for (int j = 0; j < n; j++) {
			temp2.push_back(temp1[i][j]);
		}
		temp.push_back(temp2);
	}

}

double determinant(vector<vector<double>> A, int n)
{
	double D = 0; // Initialize result 

	//  Base case : if matrix contains single element 
	if (n == 1)
		return A[0][0];

	// To store cofactors 

	int sign = 1;  // To store sign multiplier 

	 // Iterate for each element of first row 
	for (int f = 0; f < n; f++)
	{
		vector<vector<double>> temp;
		// Getting Cofactor of A[0][f] 
		getCofactor(A, temp, 0, f, n);
		D += sign * A[0][f] * determinant(temp, n - 1);

		// terms are to be added with alternate sign 
		sign = -sign;
	}

	return D;
}

// Function to get adjoint of A[N][N] in adj[N][N]. 
double** adjoint(vector<vector<double>> A)
{

	int N = A.size();
	double** adj = new double*[N];
	for (int i = 0; i < N; i++) {
		adj[i] = new double[N];
	}
	if (N == 1)
	{
		adj[0][0] = 1;
		return adj;
	}


	// temp is used to store cofactors of A[][] 
	int sign = 1;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			vector<vector<double>> temp;
			// Get cofactor of A[i][j] 
			getCofactor(A, temp, i, j, N);

			// sign of adj[j][i] positive if sum of row 
			// and column indexes is even. 
			sign = ((i + j) % 2 == 0) ? 1 : -1;

			// Interchanging rows and columns to get the 
			// transpose of the cofactor matrix 
			adj[j][i] = (sign)*(determinant(temp, N - 1));
		}
	}
	return adj;
}

// Function to calculate and store inverse, returns false if 
// matrix is singular 
bool inverse(vector<vector<double>> A, vector<vector<double>> &inverse)
{
	// Find determinant of A[][] 
	int N = A.size();
	double det = determinant(A, N);
	if (det <= 0.0001 && det >= -0.0001)
	{
		cout << "Singular matrix, can't find its inverse";
		return false;
	}

	// Find adjoint 
	double** adj = new double*[N];
	for (int i = 0; i < N; i++) {
		adj[i] = new double[N];
	}
	adj = adjoint(A);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
	for (int i = 0; i < N; i++) {
		vector<double> temp2;
		for (int j = 0; j < N; j++) {
			temp2.push_back(adj[i][j] / det);
		}
		inverse.push_back(temp2);
	}
	return true;
}


vector < vector<double> > getMatrixInverse(const vector<vector<double>>& matrix) {
	vector< vector<double> > matrixInverse;
	if (!inverse(matrix, matrixInverse))
	{
		exit(0);
	}
	return matrixInverse;
}

vector<vector<double>> getMatrixTranspose(const vector<vector<double>> &matrix) {
	vector<vector<double>> matrixTranspose;
	int m = matrix.size();
	int n = matrix[0].size();
	for (int i = 0; i < n; i++) {
		vector<double> temp;
		for (int j = 0; j < m; j++) {
			temp.push_back(matrix[j][i]);
		}
		matrixTranspose.push_back(temp);
	}
	return matrixTranspose;
}

Ray getTransformedRay(Ray r1, const vector<vector<double>> &vMat, const vector<vector<double> > &rMat) {
	Point newSource = multiplyMatrixVector(vMat, r1.getSource());
	Vector newDirection = multiplyMatrixVector(rMat, r1.getDirection());
	Ray r(newSource, newDirection);
	return r;
}

vector<vector<double>> formTranslationMatrix(double tx, double ty, double tz) {
	return { {1, 0, 0, 0}, {0,1,0,0}, {0,0,1,0}, {tx, ty, tz, 1} };
}
