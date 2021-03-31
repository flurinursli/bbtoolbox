/*
 Purpose:
   To generate random numbers based on the TRNG4 library

 Revisions:
     Date                    Description of change
     ====                    =====================
   02/09/20                  original version
   08/03/21                  leapfrog added
*/

#include <omp.h>
//#include <iostream>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/uniform_dist.hpp>
#include <trng/normal_dist.hpp>


// How to compile: g++ -c trng.cpp -Ipath-to-TRNG-include-folder -O3 -cpp (-DDOUBLE_PREC)

extern "C"
{
  void setup_trng(int seed, long skip, int streams, int nostream);
#ifdef DOUBLE_PREC
  void uni01(double a, double b, int npts, double * x);
  void uni(double a, double b, int npts, double * x);
  void norm(double mean, double std, int npts, double * x);
#else
  void uni01(float a, float b, int npts, float * x);
  void uni(float a, float b, int npts, float * x);
  void norm(float mean, float std, int npts, float * x);
#endif
}

static trng::yarn2 * r;
#pragma omp threadprivate(r)

// -----------------------------------------------------------------------------

void setup_trng(int seed, long skip, int streams, int nostream){

  unsigned long myseed = seed;

  #pragma omp parallel
  {

    delete r;

    r = new trng::yarn2;

    r->seed(myseed);
    r->jump(skip);
    r->split(streams, nostream);

  }

}



#ifdef DOUBLE_PREC
void uni01(double a, double b, int npts, double * x)
{
  trng::uniform01_dist<double> u;
#else
void uni01(float a, float b, int npts, float * x)
{
  trng::uniform01_dist<float> u;
#endif

  for (int i = 0; i < npts; i++){
    x[i] = u(*r);
  }
}

// -----------------------------------------------------------------------------

#ifdef DOUBLE_PREC
void uni(double a, double b, int npts, double * x)
{
  trng::uniform_dist<double> u(a, b);
#else
void uni(float a, float b, int npts, float * x)
{
  trng::uniform_dist<float> u(a, b);
#endif

  for (int i = 0; i < npts; i++){
    x[i] = u(*r);
  }
}

// -----------------------------------------------------------------------------

#ifdef DOUBLE_PREC
void norm(double mean, double std, int npts, double * x)
{
  trng::normal_dist<double> u(mean, std);
#else
void norm(float mean, float std, int npts, float * x)
{
  trng::normal_dist<float> u(mean, std);
#endif

  for (int i = 0; i < npts; i++){
    x[i] = u(*r);
  }
}
