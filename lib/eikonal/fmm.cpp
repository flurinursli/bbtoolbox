#include <iostream>
#include "fast_marching_method.hpp"

#ifdef DOUBLE_PREC
 using real = double;
#else
 using real = float;
#endif


extern "C"
{
  void fmm_c(int gridsize[], real gridstep[], int nsrc, int* src, real* values, real* velocity, real* arrivals);
}

//using namespace std;
namespace fmm = thinks::fast_marching_method;

void fmm_c(int gridsize[], real gridstep[], int nsrc, int* src, real* values, real* velocity, real* arrivals){

  /*
  std::cout << "gridsize: " << gridsize[0] << " " << gridsize[1] << " " <<
  "\ngridstep: " << gridstep[0] << " " << gridstep[1] << " " <<
  "\nsrc: " << src[0] << " " << src[1] << " " << std::endl;
  */

  auto npts = gridsize[0] * gridsize[1];

  //std::vector<std::array<int32_t, 2>> source_indices { {src[0]-1, src[1]-1} };
  //std::vector<real> source_times {0.f};

  // allocate memory
  std::vector<std::array<int32_t, 2>> source_indices (nsrc);
  std::vector<real> source_times (nsrc);

  // copy sources and companion initial times, subtract one to take into account that we are calling from Fortran
  for (int i = 0; i < nsrc; i++){
    source_indices[i] = {{src[2*i]-1, src[2*i+1]-1}};
    source_times[i]   = values[i];
  }

  std::array<size_t, 2> grid_size {static_cast<size_t>(gridsize[0]), static_cast<size_t>(gridsize[1])};

  std::array<real, 2> grid_spacing {gridstep[0], gridstep[1]};

  // allocate memory
  std::vector<real> speed_buffer(npts, 0.f);

  // copy velocity into std::vector
  for (int i = 0; i < npts; i++){
    speed_buffer[i] = *velocity;
    velocity++;
  }

  auto solver = fmm::HighAccuracyVaryingSpeedEikonalSolver<real, 2>(grid_spacing, grid_size, speed_buffer);

  auto arrival_times = fmm::SignedArrivalTime(grid_size, source_indices, source_times, solver);

  const auto* x = arrival_times.data();

  // copy arrival-times from std::vector
  for (int i = 0; i < gridsize[0]*gridsize[1]; i++){
    *arrivals = *x;
    x++;
    arrivals++;
  }

}
