#ifndef _BEZIER_SURFACE_H_
#define _BEZIER_SURFACE_H_
#include <math.h>
#include "Vektor.h"

class bezier_surface{
 private: 
 public: 
  bezier_surface(); 
  ~bezier_surface();
  
  virtual double operator()(double u,double v)=0;
  
}; 

#endif /* _BEZIER_SURFACE_H_ */

  
