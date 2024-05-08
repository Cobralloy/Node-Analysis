#include <header.h>

SquareAugmentedMatrix::SquareAugmentedMatrix(int equations)
{
            initialEquations = equations;
            // Intialize an augmented matrix of 0s
            matrix = vector<vector<double>> (equations, vector<double> (equations + 1, 0));
}

void SquareAugmentedMatrix::CprintMatrix() 
{
    cout << setprecision(5);
    for(int i = 0; i < matrix.size(); i++)
    {
        for(int j = 0; j < matrix[i].size(); j++)
        {
            if (j == matrix.size()) // Add the line for the constant column of the augmented matrix
                cout << " | ";

            cout << setw(10);
            if (fabs(matrix[i][j]) > 1e-9)
                cout << matrix[i][j] << " ";
            else   
                cout << 0 << " ";
        }
        cout<< endl;
    }
    cout << endl;
}

// Function to print out the augmented matrix in a nice grid in the output file
void SquareAugmentedMatrix::FprintMatrix(ofstream &fileOut) 
{
    fileOut << setprecision(5);
    for(int i = 0; i < matrix.size(); i++)
    {
        for(int j = 0; j < matrix[i].size(); j++)
        {
            if (j == matrix.size()) // Add the line for the constant column of the augmented matrix
                fileOut << " | ";

            fileOut << setw(10);
            if (fabs(matrix[i][j]) > 1e-9)
                fileOut << matrix[i][j] << " ";
            else   
                fileOut << 0 << " ";
        }
        fileOut << endl;
    }
    fileOut << endl;
}    

void SquareAugmentedMatrix::expand1()
{
    // Adds one row and one column of 0s to the matrix (keeps constant column on the very right)
    for (int row = 0; row < matrix.size(); row++)
    {
        matrix[row].insert(matrix[row].end()-1, 0.0);
    }
    matrix.push_back(vector<double> (matrix.size()+2, 0));
}

// Row and column (lead) indices start at 0
void SquareAugmentedMatrix::reduceRow(int row, int lead) 
{
    /*  Reduce row so that the leading entry is 1 by dividing every entry in the row by the leading entry
        Assumes that all entries in rowVector before the lead(th) column are 0  */
    double factor = matrix[row][lead]; 
    if (fabs(factor) > 1e-9) // Avoid divide by 0 (federal crime)
    {
        for (int entry = lead; entry < matrix[row].size(); entry++) 
        {
            matrix[row][entry] /= factor;
        }
    }
}  

// Row indices start at 0
void SquareAugmentedMatrix::subtractRow(int referenceRow, int operatingRow, int lead) 
{
    /*  Subtract a multiple of the referenceRow from the operatingRow 
        Mutiple (double factor) is such that the leading entry of operating row is eliminated by subtraction
        Assumes that all entries before the lead(th) position in the referenceRow are 0, hence int entry = lead;
        - If the entry at the lead(th) position is 0, assume that the entire row is 0s and return
            the operatingRow (subtracting by a row of 0s does nothing) */
    if (fabs(matrix[referenceRow][lead]) > 1e-9) 
    {
        // Determine multiple of referenceRow to subtract from operatingRow
        double factor = matrix[operatingRow][lead] / matrix[referenceRow][lead]; 
        for (int entry = lead; entry < matrix[referenceRow].size(); entry++) 
        {
            matrix[operatingRow][entry] -= factor * matrix[referenceRow][entry];
        }
    }
}

void SquareAugmentedMatrix::rref(ofstream &fileOut, bool debugging)
{
    // For each row, we need to make sure the leading entry is 1 by multiplying the row by a constant
    for (int referenceRow = 0; referenceRow < matrix.size(); referenceRow++) 
    {
        /*  The leading entry of the referenceRow should be on the diagonal (the nested
            for loop below makes all entries to the left 0), so we pass referenceRow to int lead */
        reduceRow(referenceRow, referenceRow);
        if (debugging)
        {
            fileOut << "Reducing row " << referenceRow + 1 << " to a leading entry of 1" << endl;
            FprintMatrix(fileOut);
        }

        /*  Once the referenceRow has a leading entry of 1, subtract a constant multiple of referenceRow
            (constant determined by the subtractRow function) from all other rows (operatingRows).
            This results in the leading 1 of the referenceRow being the only non-zero entry in its column */
        for (int operatingRow = 0; operatingRow < matrix.size(); operatingRow++) 
        {
            if (referenceRow != operatingRow) // Can't subtract a row by itself!
            {
                /*  Pass referenceRow row number to int lead because all entries before the lead(th) column 
                    in the referenceRow should be 0 */ 
                subtractRow(referenceRow, operatingRow, referenceRow);
                if (debugging)
                {
                    fileOut << "Subtracting row " << referenceRow + 1 << " from row " << operatingRow + 1 << endl; 
                    FprintMatrix(fileOut);
                }
            }
        }
    }
}