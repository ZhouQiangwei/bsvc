#ifdef _MSC_VER
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  pragma warning(disable: 4510) // default constructor could not be generated.
#  pragma warning(disable: 4610) // can never be instantiated - user defined constructor required.
#endif


#include <iostream>
#include <iomanip>
#include <boost/math/distributions/students_t.hpp>

float tstudtest(float* arr1, float* arr2, int n, int m);
