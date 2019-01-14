// CPP Program to implement t-test. 
#include <bits/stdc++.h> 
using namespace std; 
  
// Function to find mean. 
float Mean(float arr[], int n) 
{ 
    float sum = 0; 
    for (int i = 0; i < n; i++) 
        sum = sum + arr[i]; 
    return sum / n; 
} 
  
// Function to find standard 
// deviation of given array. 
float standardDeviation(float arr[], int n) 
{ 
    float sum = 0; 
    for (int i = 0; i < n; i++) 
        sum = sum + (arr[i] - Mean(arr, n)) * 
                    (arr[i] - Mean(arr, n)); 
  
    return sqrt(sum / (n - 1)); 
} 
  
// Function to find t-test of 
// two set of statistical data. 
float tTest(float arr1[], int n, 
            float arr2[], int m) 
{ 
    float mean1 = Mean(arr1, n); 
    float mean2 = Mean(arr2, m); 
    float sd1 = standardDeviation(arr1, n); 
    float sd2 = standardDeviation(arr2, m); 
  
    // Formula to find t-test 
    // of two set of data. 
    float t_test = (mean1 - mean2) / sqrt((sd1 * sd1) 
                              / n + (sd2 * sd2) / m); 
    return t_test; 
} 
  
// Driver function. 
int main() 
{ 
    float arr1[] = { 1, 1, 0, 0, 0, 0, 0, 0 }; 
  
    // Calculate size of first array. 
    int n = sizeof(arr1) / sizeof(arr1[0]); 
    float arr2[] = { 0, 0, 0, 0, 0, 0, 0, 0 }; 
  
    // Calculate size of second array. 
    int m = sizeof(arr2) / sizeof(arr2[0]); 
  
    // Function call. 
    cout << tTest(arr1, n, arr2, m); 
    
    return 0; 
} 
