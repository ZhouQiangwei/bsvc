// Copyright John Maddock 2006.
// Copyright Paul A. Bristow 2007.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef _MSC_VER
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  pragma warning(disable: 4510) // default constructor could not be generated.
#  pragma warning(disable: 4610) // can never be instantiated - user defined constructor required.
#endif

#include <iostream>
#include <iomanip>
#include <boost/math/distributions/students_t.hpp>

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

float two_samples_t_test_equal_sd(
        double Sm1,
        double Sd1,
        unsigned Sn1,
        double Sm2,
        double Sd2,
        unsigned Sn2,
        double alpha)
{
   //
   // Sm1 = Sample Mean 1.
   // Sd1 = Sample Standard Deviation 1.
   // Sn1 = Sample Size 1.
   // Sm2 = Sample Mean 2.
   // Sd2 = Sample Standard Deviation 2.
   // Sn2 = Sample Size 2.
   // alpha = Significance Level.
   //
   // A Students t test applied to two sets of data.
   // We are testing the null hypothesis that the two
   // samples have the same mean and that any difference
   // if due to chance.
   // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
   //
   using namespace std;
   using namespace boost::math;

   // Now we can calculate and output some stats:
   //
   // Degrees of freedom:
   double v = Sn1 + Sn2 - 2;
   // Pooled variance:
   double sp = sqrt(((Sn1-1) * Sd1 * Sd1 + (Sn2-1) * Sd2 * Sd2) / v);
   // t-statistic:
   double t_stat;
   if(Sm1 == Sm2) t_stat = 0;
   else if(sp ==0) t_stat = 0;
   else t_stat  = (Sm1 - Sm2) / (sp * sqrt(1.0 / Sn1 + 1.0 / Sn2));
   //
   // Define our distribution, and get the probability:
   students_t dist(v);
   double q = cdf(complement(dist, fabs(t_stat)));
   return 2*q;
}

float two_samples_t_test_unequal_sd(
        double Sm1,
        double Sd1,
        unsigned Sn1,
        double Sm2,
        double Sd2,
        unsigned Sn2,
        double alpha)
{
   //
   // Sm1 = Sample Mean 1.
   // Sd1 = Sample Standard Deviation 1.
   // Sn1 = Sample Size 1.
   // Sm2 = Sample Mean 2.
   // Sd2 = Sample Standard Deviation 2.
   // Sn2 = Sample Size 2.
   // alpha = Significance Level.
   //
   // A Students t test applied to two sets of data.
   // We are testing the null hypothesis that the two
   // samples have the same mean and that any difference
   // if due to chance.
   // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
   //
   using namespace std;
   using namespace boost::math;

   // Print header:
   // Now we can calculate and output some stats:
   //
   // Degrees of freedom:
   double v = Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2;
   v *= v;
   double t1 = Sd1 * Sd1 / Sn1;
   t1 *= t1;
   t1 /=  (Sn1 - 1);
   double t2 = Sd2 * Sd2 / Sn2;
   t2 *= t2;
   t2 /= (Sn2 - 1);
   v /= (t1 + t2);
   // t-statistic:
   double t_stat;
   if(Sm1 == Sm2) t_stat = 0;
   else if(Sd1==0 && Sd2==0) t_stat = 0;
   else t_stat  = (Sm1 - Sm2) / sqrt(Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2);
   //
   // Define our distribution, and get the probability:
   //
   students_t dist(v);
   double q = cdf(complement(dist, fabs(t_stat)));
   return 2*q;
}

float tstudtest(float* arr1, float* arr2, int n, int m)
{
   //
   // Run tests for Car Mileage sample data
   // http://www.itl.nist.gov/div898/handbook/eda/section3/eda3531.htm
   // from the NIST website http://www.itl.nist.gov.  The data compares
   // miles per gallon of US cars with miles per gallon of Japanese cars.
   //
    //float arr1[] = { 1, 1, 1, 0 };

    // Calculate size of first array.
    //int n = sizeof(arr1) / sizeof(arr1[0]);
    //float arr2[] = { 0, 0, 0, 0 };

    // Calculate size of second array.
    //int m = sizeof(arr2) / sizeof(arr2[0]);
    float mean1 = Mean(arr1, n);
    float mean2 = Mean(arr2, m);
    float sd1 = standardDeviation(arr1, n);
    float sd2 = standardDeviation(arr2, m);
    float qval=0;
    if(sd1 == sd2)
        qval = two_samples_t_test_equal_sd(mean1, sd1, n, mean2, sd2, m, 0.05);
    else qval = two_samples_t_test_unequal_sd(mean1, sd1, n, mean2, sd2, m, 0.05);
    //std::cout << n << "\t" << m << "\t" <<  sd1 << "\t" << sd2 << "\t" << qval << "\n";
   //std::cout << qval;
   //two_samples_t_test_equal_sd(20.14458, 6.414700, 249, 30.48101, 6.107710, 79, 0.05);
   //two_samples_t_test_unequal_sd(20.14458, 6.414700, 249, 30.48101, 6.107710, 79, 0.05);

   return qval;
} 

/*
Output is

------ Build started: Project: students_t_two_samples, Configuration: Debug Win32 ------
Compiling...
students_t_two_samples.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\students_t_two_samples.exe"
_______________________________________________
Student t test for two samples (equal variances)
_______________________________________________

Number of Observations (Sample 1)                      =  249
Sample 1 Mean                                          =  20.145
Sample 1 Standard Deviation                            =  6.4147
Number of Observations (Sample 2)                      =  79
Sample 2 Mean                                          =  30.481
Sample 2 Standard Deviation                            =  6.1077
Degrees of Freedom                                     =  326
Pooled Standard Deviation                              =  326
T Statistic                                            =  -12.621
Probability that difference is due to chance           =  5.273e-030

Results for Alternative Hypothesis and alpha           =  0.0500

Alternative Hypothesis              Conclusion
Sample 1 Mean != Sample 2 Mean       NOT REJECTED
Sample 1 Mean <  Sample 2 Mean       NOT REJECTED
Sample 1 Mean >  Sample 2 Mean       REJECTED


_________________________________________________
Student t test for two samples (unequal variances)
_________________________________________________

Number of Observations (Sample 1)                      =  249
Sample 1 Mean                                          =  20.14458
Sample 1 Standard Deviation                            =  6.41470
Number of Observations (Sample 2)                      =  79
Sample 2 Mean                                          =  30.48101
Sample 2 Standard Deviation                            =  6.10771
Degrees of Freedom                                     =  136.87499
T Statistic                                            =  -12.94627
Probability that difference is due to chance           =  1.571e-025

Results for Alternative Hypothesis and alpha           =  0.0500

Alternative Hypothesis              Conclusion
Sample 1 Mean != Sample 2 Mean       NOT REJECTED
Sample 1 Mean <  Sample 2 Mean       NOT REJECTED
Sample 1 Mean >  Sample 2 Mean       REJECTED

Build Time 0:03
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\students_t_two_samples\Debug\BuildLog.htm"
students_t_two_samples - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========

*/
