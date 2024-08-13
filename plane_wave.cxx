/*
 * plane_wave.cxx
 * Jordan Walsh
 * 
 */


#include <iostream>
#include <fstream>
using namespace std;
#include <string>
#include <vector>
#include <math.h>
//#include <Eigen/Core>
#include <Eigen/Eigenvalues>

const double pi = M_PI;
const double hbarSq_2m = 3.81; // [Angstrom eV]
const double m_0 = 1.0; //electron rest mass

void writeVectorToCSV(const vector<int>& vec, const string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    for (size_t i = 0; i < vec.size(); ++i) {
        file << vec[i];
        if (i < vec.size() - 1) {
            file << ",";
        }
    }
    file.close();
    std::cout << "Data written to " << filename << std::endl;
}

void writeMatrixToCSV(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            file << matrix[i][j];
            if (j < matrix[i].size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }
    file.close();
    std::cout << "Matrix written to " << filename << std::endl;
}

void print_matrix (vector<vector<double>> &M) { // square matrices only
	int n = (int)M[0].size();
	cout << "[";
	for(int i = 0; i < n; i++){	
		for(int j = 0; j < n; j++){ 
			cout << M[i][j];			//print element
			if(j != n-1) cout << ",\t";	// (decoration*)
			else if (j == n-1 && i == n-1) cout << "]";
		}
		cout << "\n ";
	}
}

double delta (int i, int j) {
	if (i == j) return 1.0;
		else return 0.0;
}

int main(int argc, char **argv)
{
	// material parameters 
	double m_b = 0.08775; // [1/m_0]
	double m_w = 0.067;	// [1/m_0]
	
	// well dimensions
	double a = 70.0; // [Angstrom]
	double L = 210.0; // [Angstrom]
	double V_0 = 0.202; // [eV]
	
	int M = 6;
	int M_ = 2*M+1;
	vector<vector<double>> H(M_, vector<double>(M_));
	
	for(int i = 0; i<M_; i++){
		for(int j = 0; j<M_; j++){
			double G_mi = 2*pi*(i-M) / L; //[1/Angstrom]
			double G_mj = 2*pi*(j-M) / L; //[1/Angstrom]
			double A = ( hbarSq_2m ) * ( G_mi*G_mj / m_b ) + V_0;
			double s_arg = (G_mi-G_mj)*a / 2.0;
			double B;
			if (i==j) B = (a/L); //Divide by 0 error otherwise, doing this produces solutions to matrix
				else B = (a/L) * ( sin(s_arg) / (s_arg) );
			double C = ( hbarSq_2m ) * (G_mi*G_mj) * ( (1.0/m_w) - (1.0/m_b) ) - V_0;
			H[i][j] = A*delta(i,j) + B*C;
		}
	}
	print_matrix(H);

	//solve Matrix (Eigen)
	Eigen::MatrixXf M_eigen(M_, M_);
	for (int i = 0; i < M_; ++i) {
		for (int j = 0; j < M_; ++j) {
			M_eigen(i, j) = H[i][j];
		}
	}

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver(M_eigen);
	Eigen::VectorXf eigenvalues = solver.eigenvalues(); // Energies
	Eigen::MatrixXf eigenvectors = solver.eigenvectors(); // Coefficients
	vector<double> outEnergies(M_);
	vector<vector<double>> outCoefficients(M_, vector<double>(M_));
	for (int i = 0; i < M_; ++i) {
		outEnergies[i] = eigenvalues(i);
		cout << "E" << i << " : " << eigenvalues(i) << endl;	
		//cout << "V: ";
		for (int j = 0; j < M_; ++j) {
			//cout << eigenvectors(j, i) << " , ";
			outCoefficients[i][j] = eigenvectors(j, i);
		}
		cout << "\n";
	}
	print_matrix(outCoefficients);
	writeMatrixToCSV(outCoefficients, "coefficients.csv");
	
	return 0;
}

