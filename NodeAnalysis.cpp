/*
    Node Analysis Program
    Alan Hu
    15 February 2023
    ReadMe.txt for instructions
*/

#include <header.h>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip> // for setw in output
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;

struct ElementProperties
{
    string type;
    double voltage;
    double current;
    double power;
};

class NodeAnalysisMatrix: public SquareAugmentedMatrix
{
    private:
        int totalNodes = initialEquations;
        
        // Find potential difference between two nodes
        double findPD (int nodeFrom, int nodeTo)
        {
            return matrix[nodeTo][matrix.size()] - matrix[nodeFrom][matrix.size()];
        }

    public:
        using SquareAugmentedMatrix::SquareAugmentedMatrix;
        
        /*  Add uncontrolled current source
            Node indices start at 0 */
        void addUCS(int nodeFrom, int nodeTo, double current)
        {
            // elementValue leaving nodeFrom
            matrix[nodeFrom].back() -= current;
            // elementValue entering nodeTo
            matrix[nodeTo].back() += current;
        }

        /*  Add uncontrolled voltage source
            Node indices start at 0 */
        void addUVS(int nodeFrom, int nodeTo, double voltage)
        {
            expand1();
            // Entries in new column: new variable is I_source
            matrix[nodeTo][matrix.size()-1] -= 1;
            matrix[nodeFrom][matrix.size()-1] += 1;

            // Entries in new row: nodeTo - nodeFrom = voltage
            matrix[matrix.size()-1][nodeTo] += 1;
            matrix[matrix.size()-1][nodeFrom] -= 1;
            matrix[matrix.size()-1][matrix.size()] += voltage;
        }

        // Node indices start at 0
        void addResistor(int nodeFrom, int nodeTo, double resistance)
        {
            /* Resistor connected to nodeFrom and nodeTo, so we add to the diagonal for both nodes */
            matrix[nodeFrom][nodeFrom] += 1/resistance;
            matrix[nodeTo][nodeTo] += 1/resistance;

            /*  Resistor between nodeFrom and nodeTo, so we subtract symmetrically */
            matrix[nodeFrom][nodeTo] -= 1/resistance;
            matrix[nodeTo][nodeFrom] -= 1/resistance;
        }

        /*  Add voltage controlled current source
            Node indices start at 0 */
        void addVCCS(int nodeFrom, int nodeTo, int controlNodeFrom, int controlNodeTo, double multiplier)
        {
            expand1();
            // Entries in new column: new variable is I_source
            matrix[nodeFrom][matrix.size()-1] += 1;
            matrix[nodeTo][matrix.size()-1] -= 1;

            // I_source = multiplier*(controlNodeTo - controlNodeFrom)
            // multiplier*controlNodeTo - multiplier*controlNodeFrom - I_source = 0
            matrix[matrix.size()-1][controlNodeFrom] -= multiplier;
            matrix[matrix.size()-1][controlNodeTo] += multiplier;
            matrix[matrix.size()-1][matrix.size()-1] -= 1;
        }

        /*  Add voltage controlled current source
            Node indices start at 0 */
        void addVCVS(int nodeFrom, int nodeTo, int controlNodeFrom, int controlNodeTo, double multiplier)
        {
            expand1();
            // Entries in new column: new variable is I_source
            matrix[nodeFrom][matrix.size()-1] += 1;
            matrix[nodeTo][matrix.size()-1] -= 1;

            /* Entries in new row
            nodeTo - nodeFrom = multiplier * (controlNodeTo - controlNodeFrom)
            nodeTo - nodeFrom - multiplier*controlNodeTo + multiplier*controlNodeFrom = 0 */
            matrix[matrix.size()-1][nodeFrom] -= 1;
            matrix[matrix.size()-1][nodeTo] += 1;
            matrix[matrix.size()-1][controlNodeFrom] += multiplier;
            matrix[matrix.size()-1][controlNodeTo] -= multiplier;
        }

        void solveCircuit(ofstream &fileOut, bool debugging)
        {
            rref(fileOut, debugging);
            for (int node = 0; node < matrix.size(); node++) 
            {
                fileOut << noshowpos;
                double entryValue = 0;
                if (fabs(matrix[node][matrix.size()]) > 1e-9)
                    entryValue = matrix[node][matrix.size()];

                if (node < totalNodes)
                {
                    fileOut << "Node " << node+1 << " voltage: " << showpos << setw(10) << entryValue << " Volts.\n";
                }
            }
        }

        // Solve voltage, current, and power of each element in the circuit
        ElementProperties solveResistor(int nodeFrom, int nodeTo, double resistance)
        {
            ElementProperties Resistor;
            Resistor.type = "Resistor";

            Resistor.voltage = findPD(nodeFrom, nodeTo);
            Resistor.current = Resistor.voltage / resistance;
            Resistor.power = Resistor.voltage * Resistor.current;

            return Resistor;
        }

        ElementProperties solveUVS(int nodeFrom, int nodeTo, double voltage, int sourceNumber)
        {
            ElementProperties VoltageSource;
            VoltageSource.type = "UVS";

            // Check that the voltage stated in the input file is the same as the difference betweeen the solved nodes
            double voltageCheck = findPD(nodeFrom, nodeTo);
            VoltageSource.voltage = voltage;
            if (fabs(VoltageSource.voltage - voltageCheck) > 1e-9)
                cout << "Voltage calculation mismatch! " << VoltageSource.voltage << " != " << voltageCheck << endl;

            // Grabs current from appropriate I_source column
            VoltageSource.current = matrix[totalNodes+sourceNumber][matrix.size()];
            VoltageSource.power = VoltageSource.voltage * VoltageSource.current;
            
            return VoltageSource;
        }

        ElementProperties solveUCS(int nodeFrom, int nodeTo, double current)
        {
            ElementProperties CurrentSource;
            CurrentSource.type = "UCS";

            CurrentSource.voltage = findPD(nodeFrom, nodeTo);
            CurrentSource.current = current;
            CurrentSource.power = CurrentSource.voltage * CurrentSource.current;

            return CurrentSource;
        }

        ElementProperties solveVCCS(int nodeFrom, int nodeTo, int controlNodeFrom, int controlNodeTo, double multiplier, int sourceNumber)
        {
            ElementProperties CurrentSource;
            CurrentSource.type = "VCCS";

            CurrentSource.voltage = findPD(nodeFrom, nodeTo);

            // Check that the solved current matches with the multiplier from the input file
            // Grabs current from appropriate I_source column
            CurrentSource.current = matrix[totalNodes+sourceNumber][matrix.size()];
            double currentCheck = multiplier*(findPD(controlNodeFrom, controlNodeTo));
            if (fabs(CurrentSource.current - currentCheck) > 1e-9)
                cout << "Current calculation mismatch! " << CurrentSource.current << " != " << currentCheck << endl;

            CurrentSource.power = CurrentSource.voltage * CurrentSource.current;

            return CurrentSource;
        }
        
        ElementProperties solveVCVS(int nodeFrom, int nodeTo, int controlNodeFrom, int controlNodeTo, double multiplier, int sourceNumber)
        {
            ElementProperties VoltageSource;
            VoltageSource.type = "VCVS";

            // Check that the solved voltage matches with the multiplier from the input file
            VoltageSource.voltage = findPD(nodeFrom, nodeTo);
            double controlVoltage = multiplier*(findPD(controlNodeFrom, controlNodeTo));
            if (fabs(VoltageSource.voltage - controlVoltage) > 1e-9)
                cout << "Voltage calculation mismatch! " << VoltageSource.voltage << " != " << controlVoltage << endl;

            // Grabs current from appropriate I_source column
            VoltageSource.current = matrix[totalNodes+sourceNumber][matrix.size()]; 
            VoltageSource.power = VoltageSource.voltage * VoltageSource.current;

            return VoltageSource;
        }
};

int main()
{
    ifstream fileIn ("./Data/circuit.txt");
    ofstream fileOut ("./Data/results.txt");
    if (!fileIn)
    {
        cout << "The file failed to open!";
        return EXIT_FAILURE;
    }
    
    int totalNodes = 0;
    bool debugging = true; // When set to true, each step of the matrix solving process will be outputted
    fileIn >> totalNodes;
    NodeAnalysisMatrix Circuit (totalNodes);

    int nodeFrom = 0, nodeTo = 0;
    int elementType = 0; 
    /*  
        0 for uncontrolled current source
        1 for uncontrolled voltage source
        2 for resistor 
        3 for voltage controlled voltage source
        4 for voltage controlled current source
    */
    double elementValue = 0.0;
    int controlNodeFrom = 0, controlNodeTo = 0, multiplier = 0;

    // Read in data from the input file
    while (fileIn >> nodeFrom >> nodeTo >> elementType >> elementValue)
    {
        // Nodes numbers start from 1 in input file
        nodeFrom--;
        nodeTo--;
        if (elementType == 0) 
        {
            Circuit.addUCS(nodeFrom, nodeTo, elementValue);
        }
        else if (elementType == 1)
        {
            Circuit.addUVS(nodeFrom, nodeTo, elementValue);
        }
        else if (elementType == 2)
        {
            Circuit.addResistor(nodeFrom, nodeTo, elementValue);
        }
        else if (elementType == 3)
        {
            fileIn >> controlNodeFrom >> controlNodeTo;
            Circuit.addVCVS(nodeFrom, nodeTo, --controlNodeFrom, --controlNodeTo, elementValue);
        }
        else if (elementType == 4)
        {
            fileIn >> controlNodeFrom >> controlNodeTo;
            Circuit.addVCCS(nodeFrom, nodeTo, --controlNodeFrom, --controlNodeTo, elementValue);
        }
        else
        {
            fileOut << "Invalid element type! Element skipped.\n";
        }
    }
    Circuit.FprintMatrix(fileOut);
    Circuit.solveCircuit(fileOut, debugging);

    // Return to start of input file
    // Now that the circuit is solved, we can go through each element again and calculate its properties
    fileIn.clear();
    fileIn.seekg(0, ios::beg);

    fileIn >> totalNodes;
    int sources = 0; // Keep track of expansions to matrix
    fileOut << setw(12) << "Element Type" << setw(16) << "Element Value" << setw(10) << "Node(-)" 
        << setw(10) << "Node(+)" << setw(14) << "Voltage(V)" << setw(14) << "Current(A)" << setw(15) << "Power(W)\n";
    while (fileIn >> nodeFrom >> nodeTo >> elementType >> elementValue)
    {
        nodeFrom--;
        nodeTo--;
        ElementProperties Element;
        switch (elementType)
        {
        case 0:
            // Uncontrolled current source doesn't expand matrix, so no incrementation to sources
            Element = Circuit.solveUCS(nodeFrom, nodeTo, elementValue);
            break;
        case 1:
            Element = Circuit.solveUVS(nodeFrom, nodeTo, elementValue, sources++);
            break;
        case 2:
            Element = Circuit.solveResistor(nodeFrom, nodeTo, elementValue);
            break;
        case 3:
            fileIn >> controlNodeFrom >> controlNodeTo;
            Element = Circuit.solveVCVS(nodeFrom, nodeTo, --controlNodeFrom, --controlNodeTo, elementValue, sources++);
            break;
        case 4:
            fileIn >> controlNodeFrom >> controlNodeTo;
            Element = Circuit.solveVCCS(nodeFrom, nodeTo, --controlNodeFrom, --controlNodeTo, elementValue, sources++);
            break;
        default:
            fileOut << "Invalid element type! Element skipped.\n";
            break;
        }
        fileOut << setw(12) << Element.type << setw(16) << elementValue << setw(10) << ++nodeFrom << setw(10) 
            << ++nodeTo << setw(14) << Element.voltage << setw(14) << Element.current << setw(14) << Element.power << endl;
    }

    fileIn.close();
    fileOut.close();
    return EXIT_SUCCESS;
}

/*
    References
    Vectors
    https://cplusplus.com/reference/vector/vector/

    Access specifiers
    https://stackoverflow.com/questions/860339/what-is-the-difference-between-public-private-and-protected-inheritance-in-c

    Inheriting constructors
    https://stackoverflow.com/questions/347358/inheriting-constructors 

    Moving the file pointer (seekg)
    https://cplusplus.com/forum/beginner/30644/ 

    Struct
    https://www.w3schools.com/cpp/cpp_structs.asp
    
*/