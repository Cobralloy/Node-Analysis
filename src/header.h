#include <cmath>
#include <string>
#include <vector>
#include <iomanip> // for setw in output
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;

class SquareAugmentedMatrix
{
    protected:
        int initialEquations = 0;
        vector<vector<double>> matrix;

    public:
        SquareAugmentedMatrix(int equations);

        // Function to print out the augmented matrix in a nice grid in the console
        void CprintMatrix() ;

        // Function to print out the augmented matrix in a nice grid in the output file
        void FprintMatrix(ofstream &fileOut) ;   

        void expand1();
        
        // Row and column (lead) indices start at 0
        void reduceRow(int row, int lead);

        // Row indices start at 0
        void subtractRow(int referenceRow, int operatingRow, int lead);
        
        void rref(ofstream &fileOut, bool debugging);
};
