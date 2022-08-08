// TOPIC: Application of Graph Theory in Cryptography

//------------------------------------------------------------
// Implemented base paper (journal):
// Name of journal: Encryption algorithm using Graph Theory
// Author: Wael Etaiwi
// Date of Publication: 7th Aug 2014
//------------------------------------------------------------

// This encryption algorithm has been implemented in C++ language

//------------------------------------------------------------
// This algorithm takes a string as input from the user.

// This input string is then passed to Encode Class which encodes the entire string according to the algorithm
// The various details during the encoding process are displayed to the user as output in the terminal
// The encoded cipher is then displayed and then shard with the Decode class
// The Decode class takes the cipher text and decodes it.

// Final output, i.e. the decoded original message is shown to the user after decoding.
//------------------------------------------------------------



// ===========================================================

// Including header files
#include<bits/stdc++.h>

// Declaring namespace
using namespace std;


// Function which multiples 2 input matrices and then returns the output matrix.
vector<vector<double>> matrixmultiplication(vector<vector<double>> m1, vector<vector<double>> m2, int n)
{
    // This 2d vector matrix stores the result of multiplication.
	vector<vector<double>> result(n,vector<double> (n,0));

    // Multiplication logic
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			for(int k=0;k<n;k++)
			{
				result[i][j]+= m1[i][k] * m2[k][j];
			}
		}
	}

    // Return the resultant matrix.
	return result;
}


// Function to print a matrix.
// This function accepts a double type 2d vector matrix and prints it.
void printMatrix(vector<vector<double>> matrix, int n)
{
    int i, j;
    for(int i=0;i<10*n+6;i++)
    {
        cout<<"+";
    }
    cout<<endl;

    // Printing logic
    for(i=0; i<n; i++)
    {
        cout<<"|";
        for(j=0; j<n; j++)
        {

            cout<<std::fixed;
            cout<<setprecision(2);
            cout<<setw(10)<<matrix[i][j];
        }
        cout<<"\t |"<<endl;
    }

    for(int i=0;i<10*n+6;i++)
    {
        cout<<"+";
    }
    cout<<endl<<endl;
}


// Class Encode
// This class contains the entire logic to encode an input string.
class Encode{

public:
    // Declaring member variables
    string message;
    string finalword;
    int N;

    // Declaring member 2d matrices
    vector<vector<double>> m1;
    vector<vector<double>> m2;
    vector<vector<double>> m3;
    vector<vector<double>> key;
    vector<vector<double>> cipher;

    // Constructor of this class
    Encode()
    {

        // Taking the input
        cout<<"Enter the message to be encrypted: ";
        cin>>message;

        // Concatenating 'A' to the beginning of the input string.
        finalword="A"+message;

        // Length of concatenated string
        N = message.length() + 1;

        // Displaying the given message.
        cout<<"The message is: "<<message<<endl;

        int i, j;

        // Initialising all 2d matrices with 0 and declaring their dimensions.
        for(i=0; i<N; i++)
        {
            vector<double> newRow (N, 0);

            m1.push_back(newRow);
            m2.push_back(newRow);
            m3.push_back(newRow);
            key.push_back(newRow);
            cipher.push_back(newRow);
        }

        // Constructing the KEY which will be used in encryption process.
        for(i=0; i<N; i++)
        {
            for(j=0; j<N; j++)
                if (i <= j)
                    key[i][j] = 1;
        }

        // Displaying the KEY.
        cout<<"\nKey generated for encryption is:"<<endl;
        printMatrix(key,N);

    }


    // Function to create the initial wieghted graph and store the weights in the adjacency matrix m1.
    void createGraph()
    {
        int i, j;

        // Storing the weights of the graph in the adjacency matrix- m1.
        for(i=1; i<N; i++)
        {
            m1[i][i-1] = double(finalword[i] - finalword[i-1]+1);
            m1[i-1][i] = m1[i][i-1];
        }

        m1[1][N-1] = double(finalword[1] - finalword[N-1]+1);
        m1[N-1][1]=m1[1][N-1];

        // Creating a complete graph
        // The constructed adjacency matrix will be converted to a complete graph by joining all remaining nodes and then adding weights (starting from 256) to all the newly constructed nodes.
        int weight = 256;

        for(i=1; i<N; i++)
        {
            for(j=1; j<N; j++)
            {
                if(i != j && m1[i][j] == 0)
                {
                    m1[i][j] = weight;
                    m1[j][i] = weight;
                    weight++;
                }
            }
        }

        // Displaying the complete graph in formof an adjacency matrix- m1
        cout<<"\nAdjacency matrix of the complete graph created is (M1):"<<endl;
        printMatrix(m1,N);
    }


    //------------------------------------------------------------
    // Functions to find the Minimum Spanning Tree of the created complete graph (m1) using PRIM's ALGORITHM

    // A utility function to find the vertex with minimum key value, from the set of vertices not yet included in MST
    int minKey(int key[], bool mstSet[])
    {
        // Initialize min value
        int min = INT_MAX, min_index;

        for (int v = 0; v < N; v++)
        {
            if (mstSet[v] == false && key[v] < min)
            {
                min = key[v], min_index = v;
            }
        }
        return min_index;
    }

    // This function converts the computed MST into an adjacency matrix
    // The adjacency matrix M2 stores the MST of the adjacency matrix M1.
    void convert_MST(int parent[])
    {
        int i,j;

        // Conversion of MST into adjacency matrix M2
        for ( i = 1; i < N; i++)
        {
            m2[i][parent[i]] = m1[i][parent[i]];
            m2[parent[i]][i]=m2[i][parent[i]];
        }

        for(i=0;i<N;i++)
        {
            m2[i][i] = i;
        }

        cout<<"\nMST formed by the above complete graph is (M2):"<<endl;


        // Displaying adjacency matrix M2 containing the MST.
        printMatrix(m2,N);


    }

    //Driver function to calculate the MST of the adjacency matrix M1
    void computeMST()
    {
        int parent[N];  // Stores the constructed MST
        int key[N];     // Key values used to pick minimum weight edge
        bool mstSet[N]; // To represent set of vertices included in MST

        // Initialising the Keys to infinite. (or max value)
        for (int i = 0; i < N; i++)
        {
            key[i] = INT_MAX, mstSet[i] = false;
        }

        // Including first 1st vertex in MST.
        // Make key 0 so that this vertex is picked as first vertex.
        key[0] = 0;
        parent[0] = -1;

        // Prim's algo main logic
        for (int count = 0; count < N - 1; count++)
        {
            // Pick the minimum key vertex from the set of vertices not yet included in MST
            int u = minKey(key, mstSet);

            // This picked key is added to mstSet
            mstSet[u] = true;

            // Updating the parent index and key value of adjacent vertices of the picked values.
            // Only those vertices are considered which have not been included yet in the MST
            for (int v = 0; v < N; v++)
            {
                if (m1[u][v] && mstSet[v] == false && m1[u][v] < key[v])
                {
                    parent[v] = u, key[v] = m1[u][v];
                }
            }
        }

        // The created MST is converted into an adjacency matrix (M1)
        convert_MST(parent);
    }

    //------------------------------------------------------------

    // This is the driver function which runs the entire Encoding process
    void superiorFunc()
    {
        cout<<"\nMessage is getting encrypted...Please wait!!\n";

        // Encoding starts...

        // Create the initial graph- adjacency matrix m1
        createGraph();

        // Compute the MST of m1 and store in matrix m2
        computeMST();

        // m3 = m1*m2
        m3 = matrixmultiplication(m1,m2,N);

        // Print m3
        cout<<"\nMarix multipliction of M1 and M2 (=M3) is:"<<endl;
        printMatrix(m3, N);

        // Calculating the cipher matrix to be sent through matrix multiplication: cipher = key * m3
        cipher = matrixmultiplication(key, m3, N);

        // Printing the matrix
        cout<<"\nCipher generated is:"<<endl;
        printMatrix(cipher, N);

        // Printing encoded message
        cout<<"\n The encoded message is : ";
        for(int i=0; i<N ; i++)
        {
            for(int j=0; j<N ; j++)
            {
                cout<<cipher[i][j]<<" ";
            }
        }
        cout<<endl;
    }


};



//Object of this class will help to decrypt the encoded message
class Decode{

public:

    // Initialising member matrices and variables
    vector<vector<double>> m1;
    vector<vector<double>> m2;
    vector<vector<double>> m3;
    vector<vector<double>> key_inverse;
    vector<vector<double>> cipher;
    int N;

    //Constructor to set the received matrices
    // The cipher matrix and m1 matrix are sent to receiver class
    Decode(int N, vector<vector<double>> m1,vector<vector<double>> cipher,vector<vector<double>> key)
    {
        // Initialising the member variables
        this->m1 = m1;
        this->cipher = cipher;
        this->key_inverse=key;
        this->N = N;

        // Displaying all the necessary matrices
        cout<<"\nMessage successfully received on receiver end...\n";
        cout<<"\nReceived Cipher is:\n"<<endl;
        printMatrix(cipher, N);
        cout<<"\nReceived graph M1 is as:\n";
        printMatrix(m1,N);
    }

    //Function to calculate determinant of a square matrix
    double getDeterminant(vector<vector<double>> vect)
    {
        if(vect.size() != vect[0].size())
        {
            throw std::runtime_error("Matrix is not quadratic");
        }
        int dimension = vect.size();

        if(dimension == 0)
        {
            return 1;
        }

        if(dimension == 1)
        {
            return vect[0][0];
        }

        //Formula for 2x2-matrix
        if(dimension == 2)
        {
            return vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
        }

        double result = 0;
        int sign = 1;
        for(int i = 0; i < dimension; i++)
        {

            //Submatrix
            vector<vector<double>> subVect(dimension - 1, vector<double> (dimension - 1));
            for(int m = 1; m < dimension; m++)
            {
                int z = 0;
                for(int n = 0; n < dimension; n++)
                {
                    if(n != i) {
                        subVect[m-1][z] = vect[m][n];
                        z++;
                    }
                }
            }

            //Recursive call
            result = result + sign * vect[0][i] * getDeterminant(subVect);
            sign = -sign;
        }

        return result;
    }

    //Function to calculate transpose of a matrix
    vector<vector<double>> getTranspose(vector<vector<double>> matrix1)
    {

        //Transpose-matrix: height = width(matrix), width = height(matrix)
        vector<vector<double>> solution(matrix1[0].size(),vector<double> (matrix1.size()));

        //Filling solution-matrix
        for(int i = 0; i < matrix1.size(); i++)
        {
            for(int j = 0; j < matrix1[0].size(); j++)
            {
                if(matrix1[i][j]==-0)
                {matrix1[i][j]=0;}
                solution[j][i] = matrix1[i][j];
            }
        }

        // Returning the result.
        return solution;
    }

    //Function to get co factor
    vector<vector<double>> getCofactor(vector<vector<double>> vect)
    {
        if(vect.size() != vect[0].size())
        {
            throw std::runtime_error("Matrix is not quadratic");
        }

        vector<vector<double>> solution(vect.size(), vector<double> (vect.size()));
        vector<vector<double>> subVect(vect.size() - 1,vector<double> (vect.size() - 1));

        for(int i = 0; i < vect.size(); i++)
        {
            for(int j = 0; j < vect[0].size(); j++)
            {

                int p = 0;
                for(int x = 0; x < vect.size(); x++)
                {
                    if(x == i)
                    {
                        continue;
                    }
                    int q = 0;

                    for(int y = 0; y < vect.size(); y++)
                    {
                        if(y == j)
                        {
                            continue;
                        }

                        subVect[p][q] = vect[x][y];
                        q++;
                    }
                    p++;
                }
                solution[i][j] = pow(-1, i + j) * getDeterminant(subVect);
            }
        }
        return solution;
    }

    //Function to get the inverse of a matrix
    vector<vector<double>> getInverse(const vector<vector<double>> vect)
    {
        // Input matrix should not be a null matrix or zero matrix.
        if(getDeterminant(vect) == 0)
        {
            throw std::runtime_error("Determinant is 0");
        }

        // Compte determinant and then its inverse.
        double d = 1.0/getDeterminant(vect);
        vector<vector<double>> solution(vect.size(), vector<double> (vect.size()));

        for(int  i = 0; i < vect.size(); i++)
        {
            for(int  j = 0; j < vect.size(); j++)
            {
                solution[i][j] = vect[i][j];
            }
        }

        // Compute Transpose
        solution = getTranspose(getCofactor(solution));

        for(int i = 0; i < vect.size(); i++)
        {
            for(int j = 0; j < vect.size(); j++)
            {
                solution[i][j] *= d;
            }
        }

        // Return inverse.
        return solution;
    }

    // This function reconstructs the original message
    void createMessage()
    {
        string message = "";
        int i, j;

        // This matrix stores the rounded m2 matrix (MST matrix)
        vector<vector<int>> newMatrix(N,vector<int>(N,0));

        // Rounding procedure- applied to matrix m2
        cout<<"\nMatrix M2 after rounding off is: "<<endl;
        for(i=0;i<10*N+4;i++)
        {
            cout<<"+";
        }
        cout<<endl;
        for(i=0; i<N; i++)
        {
            cout<<"|";
            for(j=0; j<N; j++)
            {
                double val = m2[i][j];

                // Rounding according to the sign of the value
                if(val > 0)
                {
                    newMatrix[i][j] = int(val + 0.5);
                }
                else if(val < 0)
                {
                    newMatrix[i][j] = int(val - 0.5);
                }
                else
                    newMatrix[i][j] = 0;

                cout<<setw(10)<<newMatrix[i][j];
            }
            cout<<"\t |"<<endl;
        }
        for(i=0;i<10*N+4;i++)
        {
            cout<<"+";
        }

        // Reconstructing
        i = 0;
        j = 1;

        int start = 'A';

        // Reconstructing the message
        while(newMatrix[i][j] != 0 && i<N && j<N)
        {
            start += newMatrix[i][j]-1;
            message.push_back(char(start));
            i++;
            j++;
        }

        // If the MST is not exactly cyclic
        string newString = "";
        if(i<N && j<N && newMatrix[i][j] == 0)
        {

            int start = message[0];

            start -= (newMatrix[1][N-1]-1);

            newString.push_back(char(start));

            i = N-2;
            j = N-1;

            while(newMatrix[i][j] != 0 && i>=0 && j>=0)
            {
                start -= (newMatrix[i][j]-1);
                newString.push_back(char(start));
                i--;
                j--;
            }
        }

        reverse(newString.begin(), newString.end());

        // Final message.
        message = message + newString;

        cout<<"\n\nThe decrypted message on the receiver end is: "<<message<<endl;
    }

    // This function will execute all the required functions to decrypt the encoded message
    void superiorFunc2()
    {
        cout<<"\nKey before matrix inversion is:"<<endl;
        printMatrix(key_inverse, N);

        cout<<"\nInverse of a matrix is getting fetched...Please wait...\n";
        key_inverse = getInverse(key_inverse);

        cout<<"\nKey after taking inverse is:"<<endl;
        printMatrix(key_inverse, N);

        // Reconstructing matrix m3: m3 = key_inverse * cipher
        m3 = matrixmultiplication(key_inverse, cipher, N);
        cout<<"\nMatrix M3 as decoded on receiver side is:\n"<<endl;
        printMatrix(m3, N);

        cout<<"\nInverse of a matrix is getting fetched...Please wait...\n";
        vector<vector<double>> m1_inverse = getInverse(m1);
        cout<<"\nM1 Inverse:"<<endl;
        printMatrix(m1_inverse, N);

        // Reconstructing matrix m2: m2 = m1_inverse * m3
        m2 = matrixmultiplication(m1_inverse, m3, N);
        cout<<"\nMatrix M2 as decoded on receiver side is:\n";
        printMatrix(m2, N);

        // Reconstructing the final message and printing it.
        createMessage();

        cout<<"\nProgram executed successfully!!\n";
        cout<<"\n***********************************************************************\n                               THANK YOU\n***********************************************************************\n";
    }

};


// Main Driver class
int main()
{
    cout<<"\n***********************************************************************\n                APPLICATION OF GRAPH THEORY IN CRYPTOGRAPHY\n***********************************************************************\n";

    cout<<"\nDefinitions of matrices used in the program:\n1)M1: Adjacency Matrix of the complete graph generated by taking characters in message as nodes\n2)M2: Minimum spanning tree formed by using th graph M1\n3) M3: Matrix obtained by multiplication of M1 and M2\n4)Key: It is a pre-defined shared key\n5)Cipher: Matrix obtained by multiplying marices M3 and key\n\n";

    // Creating object of class Encode and taking the message as input
    Encode ob1;
    // Encoding message
    ob1.superiorFunc();

    cout<<"\nMessage has been successfully encrypted and has been sent to the receiver side!!\n";

    // Encoded message is passed to the Decode class after creating a object of the Decode class.
    Decode ob2(ob1.N, ob1.m1, ob1.cipher,ob1.key);
    // Decoding message
    ob2.superiorFunc2();

    // End of encoding- decoding process.
}
