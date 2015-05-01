/*************************************************************************************************************
/   Pearson correlation of large-scale neuronal ensemble GCaMP Ca+2 activity
/   CorrCpp main.m
/
/   INPUT:
/   Tab-delimited .txt file with dF/F time-series
/     /Users/Max/matlab/Calcium_Imaging/FluoroSNNAP15.03.15/FluoroSNNAP_code/PearsonCpp/dF_cell.txt
/   OUTPUT:
/   Tab-delimited .txt file with Pearson correlation matrix
/     /Users/Max/matlab/Calcium_Imaging/FluoroSNNAP15.03.15/FluoroSNNAP_code/PearsonCpp/corrmatrix.txt
/
/   Tab-delimited .txt file with p-value matrix
/     /Users/Max/matlab/Calcium_Imaging/FluoroSNNAP15.03.15/FluoroSNNAP_code/PearsonCpp/pval.txt
/
/   Maximiliano Suster, 2015
**************************************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm> // for_each vector
#include <cstdlib>
#include <cmath>

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

/****************************************************************************************************
*	compares two vector lengths
*****************************************************************************************************/

template<class T>
void comp(vector<T> v1, vector<T> v2){
    if(v1.size() != v2.size()){
        cout << "vectors not the same size\n";
        exit(1);
    }
}

/*****************************************************************************************************
 *	Euclidean distance between two vectors of type T such that T has binary +,-,*
 *****************************************************************************************************/

template<class T>
double euclidean(vector<T> v1, vector<T> v2){
    comp(v1, v2);
    T diff, sum;
    
    diff = v1[0] - v2[0];
    sum = diff * diff;
    
    for (unsigned int i=1; i < v1.size(); i++){
        diff = v1[i] - v2[i];
        sum += diff * diff;
    }
    return sqrt(double(sum));
}

/*****************************************************************************************************
 *	Jaccard Coefficient.	Use for asymetric binary values
 *****************************************************************************************************/

template<class T>
double jaccard(vector<T> v1, vector<T> v2){
    comp(v1, v2);
    int f11 = 0, f00 = 0;
    
    for (unsigned int i=0; i < v1.size(); i++){
        if(v1[i] == v2[i]){
            if(v1[i])
                f11++;
            else
                f00++;
        }
    }
    return double(f11) / double(v1.size() - (f11+f00));
}

/*****************************************************************************************************
 *	Cosine
 *****************************************************************************************************/

template<class T>
double cosine(vector<T> v1, vector<T> v2){
    comp(v1, v2);
    
    T lv1 = v1[0] * v1[0];
    T lv2 = v2[0] * v2[0];
    T dot12 = v1[0] * v2[0];
    
    for (unsigned int i=0; i < v1.size(); i++){
        lv1 += v1[i] * v1[i];
        lv2 += v2[i] * v2[i];
        dot12 += v1[i] * v2[i];
    }
    return double(dot12) / ( sqrt(double(lv1)) * sqrt(double(lv2)) );
}

/*****************************************************************************************************
 *	The mean of a vector
 *****************************************************************************************************/

template<class T>
double mean(vector<T> v1){
    T sum = v1[0];
    for (unsigned int i=1; i < v1.size(); i++)
        sum += v1[i];
    return double(sum) / double(v1.size());
}

/*****************************************************************************************************
 *	The Covariance
 *****************************************************************************************************/

template<class T>
double covariance(vector<T> v1, vector<T> v2){
    comp(v1, v2);
    double mean1 = mean(v1), mean2 = mean(v2);
    double sum = (double(v1[0]) - mean1) * (double(v2[0]) - mean2);
    for (unsigned int i=1; i < v1.size(); i++){
        sum += (double(v1[i]) - mean1) * (double(v2[i]) - mean2);
    }
    return double(sum) / double(v1.size()-1);
}

/*****************************************************************************************************
 *	standard deviation the covariance where both vectors are the same.
 *****************************************************************************************************/

template<class T>
double std_dev(vector<T> v1){
    return sqrt(covariance(v1, v1));
}

/*****************************************************************************************************
 *	Pearson Correlation
 *****************************************************************************************************/

template<class T>
double pearson(vector<T> v1, vector<T> v2){
    if (std_dev(v1) * std_dev(v2) == 0){
        cout << "(a standard deviaton was 0)";
        return -2; // I dont know what to do here???
    }
    return covariance(v1,v2)/(std_dev(v1) * std_dev(v2));
}


/*****************************************************************************************************
 * compute t-statistic from Pearson r correlation
 ******************************************************************************************************/
/*
  Student´s t statistical significance can be obtained from Pearson r coefficient using:
           r x sqrt(n-2)
     t =  --------------
             sqrt(1-r)
 
  Pearson correlation is the inverse
                t
      r = ---------------
          sqrt(t^2 + n -2)
*/

double compute_ts(double r, int n)
{
    double ts = (r*sqrt(n-2))/sqrt(1-std::pow(r, 2));
    return ts;
}


/*****************************************************************************************************
 * compute p-value from t-statistic
 *****************************************************************************************************/

// df = n– 2

double compute_pvalue(double t, int df)
{
    double a = 0.36338023;
    double w = atan(t / sqrt(df));
    double s = sin(w);
    double c = cos(w);
    
    double t1, t2;
    int j1, j2, k2;
    
    if (df % 2 == 0)       // even
    {
        t1 = s;
        if (df == 2)   // special case df=2
            return (0.5 * (1 + t1));
        t2 = s;
        j1 = -1;
        j2 = 0;
        k2 = (df - 2) / 2;
    }
    else
    {
        t1 = w;
        if (df == 1)            // special case df=1
            return 1 - (0.5 * (1 + (t1 * (1 - a))));
        t2 = s * c;
        t1 = t1 + t2;
        if (df == 3)            // special case df=3
            return 1 - (0.5 * (1 + (t1 * (1 - a))));
        j1 = 0;
        j2 = 1;
        k2 = (df - 3)/2;
    }
    for (int i=1; i>= k2; i++)
    {
        j1 = j1 + 2;
        j2 = j2 + 2;
        t2 = t2 * c * c * j1/j2;
        t1 = t1 + t2;
    }
    return 1 - (0.5 * (1 + (t1 * (1 - a * (df % 2)))));
}


// declare read data function
int read_data();


// create a new type: matrix of std::vector<float>
typedef std::vector<std::vector<float>> matrix;

/*****************************************************************************************************
* intialization of global variables
*****************************************************************************************************/

matrix dF_cell;              // matrix containing the dF/F values
std::vector<float> vec;      // a vector for storing each cell dF/F time series
int counter = 0;             // a counter for keeping track of the time-series


// source directory and files
const string rootfolder = "/Users/Max/Matlab/Calcium_Imaging/FluoroSNNAP15.03.15/FluoroSNNAP_code/PearsonCpp/";
const string infile = rootfolder + "dF_cell.txt";
const string outfile1 = rootfolder + "corrmatrix.txt";
const string outfile2 = rootfolder + "pval.txt";


/*****************************************************************************************************
 * read data function
 *****************************************************************************************************/

int read_data()
{
    
    try {
        
        std::ifstream myfile(infile);
        float val=0.0;
        std::string line;
        
        counter = 0;
        
        while (std::getline(myfile, line))
        {
            counter++;
            std::istringstream iss(line);
            //val=atof(line.c_str());
            
            while(iss >> val)
            {
                vec.push_back(val);
            }
            
            dF_cell.push_back(vec);
            vec.clear();
        }
        
        myfile.close();
        return 1;
    }
    
    catch(...)
    {
        std::cout << "Invalid input file format." << std::endl;
        return -2;  // error
    }
    
}

/*****************************************************************************************************/


int main() {
    
    // reading dF/F values from text file (exported from Matlab)
    // each line is passed into an std::vector, then stored in a matrix
    read_data();
    
    // open file for writing results to a text file
    ofstream fcorrmatrix(outfile1);
    ofstream fpval(outfile2);
    
    // measure execution of computation
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    // declare n and degrees of freedom
    int n = int(dF_cell.size());
    int df =  n - 2;
    
    if (fcorrmatrix.is_open() && fpval.is_open())
    {
        counter = 0;
        for (matrix::iterator it1 = dF_cell.begin(); it1 != dF_cell.end(); ++it1) {
            counter++;
            for (matrix::iterator it2 = dF_cell.begin(); it2 != dF_cell.end(); ++it2) {
                 double r = pearson(*it1,*it2);
                 fcorrmatrix << r << "\t";
                 double ts = compute_ts(r, n);
                 double pval = compute_pvalue(ts, df);
                 fpval << pval << "\t";
            }
            fcorrmatrix << "\n";
            fpval << "\n";
        }
     
       std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
       auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
       std::cout << "Execution time: " << float(duration/1000000) << " seconds. \n";

        // close files when finished
        // fcorrmatrix.close();
        fpval.close();
        
    }
    
    else {
    
        cout << "Unable to open output file";
        return -1;   // error
    }
    
    // report how many series were used for computing the Pearson correlation
    std::cout << "Saved Pearson correlation and p-value matrices for " << counter << " GCaMP+ neurons." << std::endl;
    
    return 0;
}


