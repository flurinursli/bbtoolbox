#include <string>
#include <proj_api.h>
//#include <iostream>

extern "C"{
void fun_r(float x, float y, float lon_0, float lat_0, float* lon, float* lat);
}

void fun_r(float x, float y, float lon_0, float lat_0, float* lon, float* lat){

  projPJ pj_longlat, pj_utm;

  std::string proj = "+proj=utm +datum=WGS84 +units=m +lon_0=" + std::to_string(lon_0) + " +lat_0=" + std::to_string(lat_0);

  pj_utm = pj_init_plus(proj.c_str());

  pj_longlat = pj_init_plus("+proj=longlat +datum=WGS84");

  double x0 = x;
  double y0 = y;

  int p = pj_transform(pj_utm, pj_longlat, 1, 1, &x0, &y0, NULL);

  *lon = (float) x0 * RAD_TO_DEG;
  *lat = (float) y0 * RAD_TO_DEG;

  //std::cout << lon << ' ' << lat << ' ' << lon_0 << ' ' << lat_0 << ' ' << y0 << ' ' << y <<  std::endl;
  //std::cout << proj.c_str() << std::endl;

  pj_free(pj_utm);
  pj_free(pj_longlat);


}
