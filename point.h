#pragma once

#include <cmath>
#include <vector>
using namespace std;

class Vector{
    public:
        double i;
        double j;
        double k;
        Vector(double a, double b, double c):i(a), j(b), k(c){}
        Vector operator +(const Vector&v){
            Vector result = v;
            result.i += this->i;
            result.j += this->j;
            result.k += this->k;
            return result;
        }
};

class Point{
    public:
        double x;
        double y;
        double z;
        Point(double a, double b, double c):x(a), y(b), z(c){}
};

Vector getUnitVector(Vector v) {
	double norm = sqrt((v.i)*(v.i) + (v.j)*(v.j) + (v.k)*(v.k));
	v.i /= norm;
	v.j /= norm;
	v.k /= norm;
	return v;
}

vector<double> multiplyVectorsPointwise(vector<double>v1, vector<double>v2) {
	// Both vectors should be of same size
	for (int i = 0;i < v1.size();i++) {
		v1[i] *= v2[i];
	}
	return v1;
}

vector<double> multiplyVectorDouble(double a, vector<double> v) {
	for (int i = 0;i < v.size();i++) {
		v[i] *= a;
	}
	return v;
}

Vector multiplyVectorDouble(double a, Vector v) {
	v.i *= a;
	v.j *= a;
	v.k *= a;
	return v;
}

Point addPointVector(Point p, Vector v) {
	p.x += v.i;
	p.y += v.j;
	p.z += v.k;
	return p;
}

Vector crossProduct(Point p1, Point p2, Point p3) {
	// forms two vectors u(p2p1) and v(p3p2) and returns u cross v
	double u1 = p1.x - p2.x;
	double u2 = p1.y - p2.y;
	double u3 = p1.z - p2.z;
	double v1 = p2.x - p3.x;
	double v2 = p2.y - p3.y;
	double v3 = p2.z - p3.z;
	double uvi = u2 * v3 - v2 * u3;
	double uvj = v1 * u3 - u1 * v3;
	double uvk = u1 * v2 - v1 * u2;
	double norm = sqrt(uvi*uvi + uvj * uvj + uvk * uvk);
	Vector v(uvi / norm, uvj / norm, uvk / norm);
	return v;
}

double dotProduct(Vector p1, Vector p2) {
	double uv = p1.i*p2.i + p1.j*p2.j + p1.k*p2.k;
	return uv;
}

Vector getSubtractionVector(Point p1, Point p2) {
	// vector p1p2 (p2 - p1)
	double x = p2.x - p1.x;
	double y = p2.y - p1.y;
	double z = p2.z - p1.z;
	Vector v(x, y, z);
	return v;
}

Point multiplyMatrix(const vector<vector<double>> &matrix, Point p) {
	// multiply a matrix(4*4) with an augmented point (input in 3d form only) to get the point in the new coordinate system
	double u1 = p.x*matrix[0][0] + p.y*matrix[1][0] + p.z*matrix[2][0] + 1*matrix[3][0];
	double u2 = p.x*matrix[0][1] + p.y*matrix[1][1] + p.z*matrix[2][1] + 1*matrix[3][1];
	double u3 = p.x*matrix[0][2] + p.y*matrix[1][2] + p.z*matrix[2][2] + 1*matrix[3][2];
	Point r(u1, u2, u3);
	return r;
}

Vector multiplyMatrix(const vector<vector<double>> &matrix, Vector v) {
	// multiply a matrix(4*4) with an augmented direction (input in 3d form only) to get the point in the new coordinate system
	double u1 = v.i*matrix[0][0] + v.j*matrix[1][0] + v.k*matrix[2][0] + 0 * matrix[3][0];
	double u2 = v.i*matrix[0][1] + v.j*matrix[1][1] + v.k*matrix[2][1] + 0 * matrix[3][1];
	double u3 = v.i*matrix[0][2] + v.j*matrix[1][2] + v.k*matrix[2][2] + 0 * matrix[3][2];
	Vector d(u1, u2, u3);
	return d;
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

	vector<vector<double>> temp; // To store cofactors 

	int sign = 1;  // To store sign multiplier 

	 // Iterate for each element of first row 
	for (int f = 0; f < n; f++)
	{
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

class Ray{
    Point source;
    Vector direction;
    public:
        Ray(Point p, Vector v) : source(p), direction(getUnitVector(v)){}
        Point getSource(){
            return source;
        }
        Vector getDirection(){
            return direction;
        }
};

class IntersectionPoint{
    Point location; // location of the intersection
    Vector normal; // normal at the location of intersection
    double rayParameter; //p0 + tp1 // tha value of t
    public:
        IntersectionPoint(Point p, Vector n, double r) : location(p), normal(n), rayParameter(r){}
        Point getLocation(){
            return location;
        }
        Vector getNormal(){
            return normal;
        }
        double getRayParameter(){
            return rayParameter;
        }
};

Ray getTransformedRay(Ray r1, const vector<vector<double> > &matrix) {
	Point newSource = multiplyMatrix(matrix, r1.getSource());
	Vector newDirection = multiplyMatrix(matrix, r1.getDirection());
	Ray r(newSource, newDirection);
	return r;
}