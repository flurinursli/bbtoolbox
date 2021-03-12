#include <string>
#include <proj_api.h>
//#include <iostream>

extern "C"{
void fun_c(float lon, float lat, float lon_0, float lat_0, float* x, float* y);
}

void fun_c(float lon, float lat, float lon_0, float lat_0, float* x, float* y){

   projPJ pj_proj, pj_longlat;

   pj_longlat = pj_init_plus("+proj=longlat +datum=WGS84");

   std::string proj = "+units=m +proj=utm +datum=WGS84 +lon_0=" + std::to_string(lon_0) + " +lat_0=" + std::to_string(lat_0);

   pj_proj = pj_init_plus(proj.c_str());

   double x0 = lon * DEG_TO_RAD;
   double y0 = lat * DEG_TO_RAD;

   int p = pj_transform(pj_longlat, pj_proj, 1, 1, &x0, &y0, NULL);

   *x = (float) x0;
   *y = (float) y0;

   //std::cout << lon << ' ' << lat << ' ' << lon_0 << ' ' << lat_0 << ' ' << y0 << ' ' << y <<  std::endl;
   //std::cout << proj.c_str() << std::endl;

   pj_free(pj_longlat);
   pj_free(pj_proj);

}
