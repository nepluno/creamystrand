#ifndef ARRAY3_UTILS_H
#define ARRAY3_UTILS_H

#include "array3.hh"
#include "util.hh"
namespace strandsim{
namespace bridson
{

template<class T>
inline Eigen::Matrix<T, 3, 1> grad_trilerp(const T& v000, const T& v100,
                                         const T& v010, const T& v110,
                                         const T& v001, const T& v101,
                                         const T& v011, const T& v111,
                                         T fx, T fy, T fz)
{
  return
  Eigen::Matrix<T, 3, 1>(-(fy - 1.) * (fz - 1.), -(fx - 1.) * (fz - 1.), -(fx - 1.) * (fy - 1.)) * v000 +
  Eigen::Matrix<T, 3, 1>((fy - 1.) * (fz - 1.), fx * (fz - 1.), fx * (fy - 1.)) * v100 +
  Eigen::Matrix<T, 3, 1>(fy * (fz - 1.), (fx - 1.) * (fz - 1.), fy * (fx - 1.) ) * v010 +
  Eigen::Matrix<T, 3, 1>(-fy * (fz - 1.), -fx * (fz - 1.), -fx * fy ) * v110 +
  Eigen::Matrix<T, 3, 1>(fz * (fy - 1.), fz * (fx - 1.), (fx - 1.) * (fy - 1.)) * v001 +
  Eigen::Matrix<T, 3, 1>(-fz * (fy - 1.), -fx * fz, -fx * (fy - 1.)) * v101 +
  Eigen::Matrix<T, 3, 1>(-fy * fz, -fz * (fx - 1.), -fy * (fx - 1.) ) * v011 +
  Eigen::Matrix<T, 3, 1>(fy * fz, fx * fz, fx * fy ) * v111;
}

template<class T>
Eigen::Matrix<T, 3, 1> affine_interpolate_value(const Eigen::Matrix<T, 3, 1>& point, const Array3<T, Array1<T> >& grid) {
  int i,j,k;
  T fx,fy,fz;
  
  get_barycentric(point[0], i, fx, 0, grid.ni);
  get_barycentric(point[1], j, fy, 0, grid.nj);
  get_barycentric(point[2], k, fz, 0, grid.nk);
  
  return grad_trilerp(
                     grid(i,j,k), grid(i+1,j,k),
                     grid(i,j+1,k), grid(i+1,j+1,k),
                     grid(i,j,k+1), grid(i+1,j,k+1),
                     grid(i,j+1,k+1), grid(i+1,j+1,k+1),
                     fx, fy, fz);
}

void accumulate_value(
   const Vec3x& point, 
   const Scalar& val, 
   Scalar* grid, int ni, int nj, int nk) {

   int i,j,k;
   Scalar fi,fj,fk;

   get_barycentric(point[0], i, fi, 0, ni);
   get_barycentric(point[1], j, fj, 0, nj);
   get_barycentric(point[2], k, fk, 0, nk);   

   auto index = [&](int i, int j, int k) -> int {
      return k * ni * nj + j * ni + i;
   };

   grid[index(i, j, k)] += (val) * (1.0 - fi) * (1.0 - fj) * (1.0 - fk);
   grid[index(i + 1, j, k)] += (val) * fi * (1.0 - fj) * (1.0 - fk);   
   grid[index(i, j + 1, k)] += (val) * (1.0 - fi) * fj * (1.0 - fk);
   grid[index(i + 1, j + 1, k)] += (val) * fi * fj * (1.0 - fk);   
   grid[index(i, j, k + 1)] += (val) * (1.0 - fi) * (1.0 - fj) * fk;
   grid[index(i + 1, j, k + 1)] += (val) * fi * (1.0 - fj) * fk;   
   grid[index(i, j + 1, k + 1)] += (val) * (1.0 - fi) * fj * fk;
   grid[index(i + 1, j + 1, k + 1)] += (val) * fi * fj * fk;   
}

void accumulate_value(
   const Vec3x& point, 
   const Scalar& val, 
   const Vec3x& c,
   Array3<Scalar, Array1<Scalar> >& grid, 
   Array3<Scalar, Array1<Scalar> >& grid_w,
   const Scalar& dx) {

   int i,j,k;
   Scalar fi,fj,fk;

   get_barycentric(point[0], i, fi, 0, grid.ni);
   get_barycentric(point[1], j, fj, 0, grid.nj);
   get_barycentric(point[2], k, fk, 0, grid.nk);   

   grid(i, j, k) += (val + c.dot(point - Vec3x(i, j, k)) * dx) * (1.0 - fi) * (1.0 - fj) * (1.0 - fk);
   grid(i + 1, j, k) += (val + c.dot(point - Vec3x(i + 1, j, k)) * dx) * fi * (1.0 - fj) * (1.0 - fk);   
   grid(i, j + 1, k) += (val + c.dot(point - Vec3x(i, j + 1, k)) * dx) * (1.0 - fi) * fj * (1.0 - fk);
   grid(i + 1, j + 1, k) += (val + c.dot(point - Vec3x(i + 1, j + 1, k)) * dx) * fi * fj * (1.0 - fk);   
   grid(i, j, k + 1) += (val + c.dot(point - Vec3x(i, j, k + 1)) * dx) * (1.0 - fi) * (1.0 - fj) * fk;
   grid(i + 1, j, k + 1) += (val + c.dot(point - Vec3x(i + 1, j, k + 1)) * dx) * fi * (1.0 - fj) * fk;   
   grid(i, j + 1, k + 1) += (val + c.dot(point - Vec3x(i, j + 1, k + 1)) * dx) * (1.0 - fi) * fj * fk;
   grid(i + 1, j + 1, k + 1) += (val + c.dot(point - Vec3x(i + 1, j + 1, k + 1)) * dx) * fi * fj * fk;   

   grid_w(i, j, k) += (1.0 - fi) * (1.0 - fj) * (1.0 - fk);
   grid_w(i + 1, j, k) += fi * (1.0 - fj) * (1.0 - fk);   
   grid_w(i, j + 1, k) += (1.0 - fi) * fj * (1.0 - fk);
   grid_w(i + 1, j + 1, k) += fi * fj * (1.0 - fk);   
   grid_w(i, j, k + 1) += (1.0 - fi) * (1.0 - fj) * fk;
   grid_w(i + 1, j, k + 1) += fi * (1.0 - fj) * fk;   
   grid_w(i, j + 1, k + 1) += (1.0 - fi) * fj * fk;
   grid_w(i + 1, j + 1, k + 1) += fi * fj * fk;   
}

Scalar interpolate_value(const Vec3x& point, const Scalar* grid, int ni, int nj, int nk) {
   int i,j,k;
   Scalar fi,fj,fk;

   get_barycentric(point[0], i, fi, 0, ni);
   get_barycentric(point[1], j, fj, 0, nj);
   get_barycentric(point[2], k, fk, 0, nk);

   auto index = [&](int i, int j, int k) -> int {
      return k * ni * nj + j * ni + i;
   };

   return trilerp(
         grid[index(i,j,k)], grid[index(i+1,j,k)], grid[index(i,j+1,k)], grid[index(i+1,j+1,k)], 
         grid[index(i,j,k+1)], grid[index(i+1,j,k+1)], grid[index(i,j+1,k+1)], grid[index(i+1,j+1,k+1)], 
         fi,fj,fk);
}
    
template<class T>
T interpolate_value(const Vec3x& point, const Array3<T, Array1<T> >& grid) {
    int i,j,k;
    T fi,fj,fk;
    
    get_barycentric(point[0], i, fi, 0, grid.ni);
    get_barycentric(point[1], j, fj, 0, grid.nj);
    get_barycentric(point[2], k, fk, 0, grid.nk);
    
    return trilerp(
                   grid(i,j,k), grid(i+1,j,k), grid(i,j+1,k), grid(i+1,j+1,k),
                   grid(i,j,k+1), grid(i+1,j,k+1), grid(i,j+1,k+1), grid(i+1,j+1,k+1),
                   fi,fj,fk);
}

Scalar interpolate_gradient(Vec3x& gradient, const Vec3x& point, const Array3x& grid) {
   int i,j,k;
   Scalar fx,fy,fz;
   
   get_barycentric(point[0], i, fx, 0, grid.ni);
   get_barycentric(point[1], j, fy, 0, grid.nj);
   get_barycentric(point[2], k, fz, 0, grid.nk);
   
   Scalar v000 = grid(i,j,k);
   Scalar v001 = grid(i,j,k+1);
   Scalar v010 = grid(i,j+1,k);
   Scalar v011 = grid(i,j+1,k+1);
   Scalar v100 = grid(i+1,j,k);
   Scalar v101 = grid(i+1,j,k+1);
   Scalar v110 = grid(i+1,j+1,k);
   Scalar v111 = grid(i+1,j+1,k+1);

   Scalar ddx00 = (v100 - v000);
   Scalar ddx10 = (v110 - v010);
   Scalar ddx01 = (v101 - v001);
   Scalar ddx11 = (v111 - v011);
   Scalar dv_dx = bilerp(ddx00,ddx10,ddx01,ddx11, fy,fz);

   Scalar ddy00 = (v010 - v000);
   Scalar ddy10 = (v110 - v100);
   Scalar ddy01 = (v011 - v001);
   Scalar ddy11 = (v111 - v101);
   Scalar dv_dy = bilerp(ddy00,ddy10,ddy01,ddy11, fx,fz);

   Scalar ddz00 = (v001 - v000);
   Scalar ddz10 = (v101 - v100);
   Scalar ddz01 = (v011 - v010);
   Scalar ddz11 = (v111 - v110);
   Scalar dv_dz = bilerp(ddz00,ddz10,ddz01,ddz11, fx,fy);

   gradient[0] = dv_dx;
   gradient[1] = dv_dy;
   gradient[2] = dv_dz;
   
   //return value for good measure.
   return trilerp(
      v000, v100,
      v010, v110, 
      v001, v101,
      v011, v111,
      fx, fy, fz);
}
}
}

#endif
