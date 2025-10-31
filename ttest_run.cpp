// CPP Program to implement t-test.
//
#include <bits/stdc++.h>
#include "tstudenttest.hpp"
using namespace std;


int main()
{
    float arr1[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1} ;//{ 1, 1, 0, 0, 0, 0, 0, 0 };

    // Calculate size of first array.
    int n = sizeof(arr1) / sizeof(arr1[0]);
    float arr2[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    // Calculate size of second array.
    int m = sizeof(arr2) / sizeof(arr2[0]);

    int totaldepth = 8, ad = 0;
    float array1[totaldepth];
    for(int i=0; i<totaldepth; i++) array1[i]=0;
    for(int i=0; i<ad; i++){
        array1[i]=1;
    }
    float array2[totaldepth];
    for(int i=0; i<totaldepth; i++) array2[i]=0;

    // Function call.
    float valu =  tstudtest(array1, array2, totaldepth, totaldepth); //tTest(array1, totaldepth, array2, totaldepth); //tTest(arr1, n, arr2, m);
    cout << valu; //tTest(arr1, n, arr2, m);

    return 0;
}
