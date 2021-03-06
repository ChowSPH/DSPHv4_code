/*
 <DUALSPHYSICS>  Copyright (c) 2015, Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

//#############################################################################
//# Descripcion:
//# =============
//# Conjunto de funciones tipicas de geometria y demas.
//#
//# Cambios:
//# =========
//# - Implementacion. (19/03/2013)
//# - Metodos para calcular area de un triangulo. (01/04/2013)
//# - Nuevos metodos para interpolacion lineal y bilineal. (08/05/2013)
//# - Nuevas funciones trigonometricas. (20/08/2015)
//#############################################################################


/// \file FunctionsMath.h \brief Declares basic/general math functions.

#ifndef _FunctionsMath_
#define _FunctionsMath_

#include "TypesDef.h"
#include <cstdlib>
#include <cmath>

/// Implements a set of basic/general math functions.
namespace fmath{

//==============================================================================
// Devuelve la interpolacion lineal de dos valores.
//==============================================================================
inline double InterpolationLinear(double x,double x0,double x1,double v0,double v1){
  const double fx=(x-x0)/(x1-x0);
  return(fx*(v1-v0)+v0);
}

//==============================================================================
// Devuelve la interpolacion bilineal de cuatro valores que forman un cuadrado.
//==============================================================================
inline double InterpolationBilinear(double x,double y,double px,double py,double dx,double dy,double vxy,double vxyy,double vxxy,double vxxyy){
  double vy0=InterpolationLinear(x,px,px+dx,vxy,vxxy);
  double vy1=InterpolationLinear(x,px,px+dx,vxyy,vxxyy);
  return(InterpolationLinear(y,py,py+dy,vy0,vy1));
}


//==============================================================================
/// Devuelve el producto escalar de 2 vectores.
//==============================================================================
inline double ProductScalar(tdouble3 v1,tdouble3 v2){
  return(v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}

//==============================================================================
/// Devuelve el producto escalar de 2 vectores.
//==============================================================================
inline float ProductScalar(tfloat3 v1,tfloat3 v2){
  return(v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}


//==============================================================================
/// Devuelve el producto vectorial de 2 vectores.
//==============================================================================
inline tdouble3 ProductVec(const tdouble3 &v1,const tdouble3 &v2){
  tdouble3 r;
  r.x=v1.y*v2.z - v1.z*v2.y;
  r.y=v1.z*v2.x - v1.x*v2.z;
  r.z=v1.x*v2.y - v1.y*v2.x;
  return(r);
}

//==============================================================================
/// Devuelve el producto vectorial de 2 vectores.
//==============================================================================
inline tfloat3 ProductVec(const tfloat3 &v1,const tfloat3 &v2){
  tfloat3 r;
  r.x=v1.y*v2.z - v1.z*v2.y;
  r.y=v1.z*v2.x - v1.x*v2.z;
  r.z=v1.x*v2.y - v1.y*v2.x;
  return(r);
}


//==============================================================================
/// Devuelve la distancia entre un punto y un plano.
//==============================================================================
inline double DistPlane(const tdouble4 &pla,const tdouble3 &pt){ 
  return(fabs(pla.x*pt.x+pla.y*pt.y+pla.z*pt.z+pla.w)/sqrt(pla.x*pla.x+pla.y*pla.y+pla.z*pla.z));
}

//==============================================================================
/// Devuelve la distancia entre un punto y un plano.
//==============================================================================
inline float DistPlane(const tfloat4 &pla,const tfloat3 &pt){ 
  return(fabs(pla.x*pt.x+pla.y*pt.y+pla.z*pt.z+pla.w)/sqrt(pla.x*pla.x+pla.y*pla.y+pla.z*pla.z));
}


//==============================================================================
/// Devuelve la distancia entre dos puntos.
//==============================================================================
inline double DistPoints(const tdouble3 &p1,const tdouble3 &p2){
  const tdouble3 v=p1-p2;
  return(sqrt(v.x*v.x+v.y*v.y+v.z*v.z));
}

//==============================================================================
/// Devuelve la distancia entre dos puntos.
//==============================================================================
inline float DistPoints(const tfloat3 &p1,const tfloat3 &p2){
  const tfloat3 v=p1-p2;
  return(sqrt(v.x*v.x+v.y*v.y+v.z*v.z));
}


//==============================================================================
/// Devuelve la distancia al (0,0,0).
//==============================================================================
inline double DistPoint(const tdouble3 &p1){
  return(sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z));
}

//==============================================================================
/// Devuelve la distancia al (0,0,0).
//==============================================================================
inline float DistPoint(const tfloat3 &p1){
  return(sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z));
}


//==============================================================================
/// Devuelve vector unitario del vector (0,0,0)->p1.
//==============================================================================
inline tdouble3 VecUnitary(const tdouble3 &p1){
  return(p1/TDouble3(DistPoint(p1)));
}

//==============================================================================
/// Devuelve vector unitario del vector (0,0,0)->p1.
//==============================================================================
inline tfloat3 VecUnitary(const tfloat3 &p1){
  return(p1/TFloat3(DistPoint(p1)));
}


//==============================================================================
/// Devuelve la normal de un triangulo.
//==============================================================================
inline tdouble3 NormalTriangle(const tdouble3& p1,const tdouble3& p2,const tdouble3& p3){
  return(ProductVec(p1-p2,p2-p3));
}

//==============================================================================
/// Devuelve la normal de un triangulo.
//==============================================================================
inline tfloat3 NormalTriangle(const tfloat3& p1,const tfloat3& p2,const tfloat3& p3){
  return(ProductVec(p1-p2,p2-p3));
}


//==============================================================================
/// Calcula el determinante de una matriz de 3x3.
//==============================================================================
inline double Determinant3x3(const tmatrix3d &d){
  return(d.a11 * d.a22 * d.a33 + d.a12 * d.a23 * d.a31 + d.a13 * d.a21 * d.a32 - d.a31 * d.a22 * d.a13 - d.a32 * d.a23 * d.a11 - d.a33 * d.a21 * d.a12);
}

//==============================================================================
/// Calcula el determinante de una matriz de 3x3.
//==============================================================================
inline float Determinant3x3(const tmatrix3f &d){
  return(d.a11 * d.a22 * d.a33 + d.a12 * d.a23 * d.a31 + d.a13 * d.a21 * d.a32 - d.a31 * d.a22 * d.a13 - d.a32 * d.a23 * d.a11 - d.a33 * d.a21 * d.a12);
}


//==============================================================================
/// Devuelve el plano formado por 3 puntos.
//==============================================================================
tdouble4 Plane3Pt(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3);

//==============================================================================
/// Devuelve el plano formado por 3 puntos.
//==============================================================================
tfloat4 Plane3Pt(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3);


//==============================================================================
/// Devuelve los tres planos normales que limitan un triangulo formado por 3 puntos.
/// Con openingdist puedes abrir o cerrar los planos normales.
//==============================================================================
void NormalPlanes3Pt(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3,double openingdist,tdouble4 &pla1,tdouble4 &pla2,tdouble4 &pla3);

//==============================================================================
/// Devuelve los tres planos normales que limitan un triangulo formado por 3 puntos.
/// Con openingdist puedes abrir o cerrar los planos normales.
/// Los calculos internos se hacen con double precision.
//==============================================================================
inline void NormalPlanes3Pt_dbl(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3,float openingdist,tfloat4 &pla1,tfloat4 &pla2,tfloat4 &pla3){
  tdouble4 plad1,plad2,plad3;
  NormalPlanes3Pt(ToTDouble3(p1),ToTDouble3(p2),ToTDouble3(p3),double(openingdist),plad1,plad2,plad3);
  pla1=ToTFloat4(plad1); pla2=ToTFloat4(plad2); pla3=ToTFloat4(plad3);
}

//==============================================================================
/// Devuelve los tres planos normales que limitan un triangulo formado por 3 puntos.
/// Con openingdist puedes abrir o cerrar los planos normales.
//==============================================================================
void NormalPlanes3Pt(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3,float openingdist,tfloat4 &pla1,tfloat4 &pla2,tfloat4 &pla3);


//==============================================================================
/// Devuelve punto de interseccion entre 3 planos no paralelos entre si.
//==============================================================================
tdouble3 Intersec3Planes(const tdouble4 &pla1,const tdouble4 &pla2,const tdouble4 &pla3);

//==============================================================================
/// Devuelve punto de interseccion entre 3 planos no paralelos entre si.
//==============================================================================
tfloat3 Intersec3Planes(const tfloat4 &pla1,const tfloat4 &pla2,const tfloat4 &pla3);


//==============================================================================
/// A partir de un triangulo formado por 3 puntos devuelve los puntos que forman
/// un triangulo mas o menos abierto segun openingdist.
//==============================================================================
void OpenTriangle3Pt(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3,double openingdist,tdouble3 &pt1,tdouble3 &pt2,tdouble3 &pt3);

  //==============================================================================
/// A partir de un triangulo formado por 3 puntos devuelve los puntos que forman
/// un triangulo mas o menos abierto segun openingdist.
//==============================================================================
void OpenTriangle3Pt(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3,float openingdist,tfloat3 &pt1,tfloat3 &pt2,tfloat3 &pt3);

//==============================================================================
/// Devuelve el area de un triangulo formado por 3 puntos.
//==============================================================================
double AreaTriangle(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3);

//==============================================================================
/// Devuelve el area de un triangulo formado por 3 puntos.
//==============================================================================
float AreaTriangle(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3);


////==============================================================================
///// Devuelve el angulo en grados que forman dos vectores.
////==============================================================================
//double JSpaceDraw::AngleVector(tdouble3 v1,tdouble3 v2){
//  return(acos(ProductScalar(v1,v2)/(ModuloVector(v1)*ModuloVector(v2)))*TODEG);
//}
////==============================================================================
///// Devuelve el plano formado por 1 vector y un punto.
////==============================================================================
//jplane JSpaceDraw::PlaneVect1Pt(tdouble3 v,tdouble3 p){ 
//  jplane plano;
//  plano.a=v.x; plano.b=v.y; plano.c=v.z; 
//  plano.d=-((v.x*p.x)+(v.y*p.y)+(v.z*p.z));
//  return(plano);
//}



//==============================================================================
/// Returns cotangent of angle in radians.
//==============================================================================
inline double cot(double z){ return(1.0 / tan(z)); }

//==============================================================================
/// Returns hyperbolic cotangent of angle in radians.
//==============================================================================
inline double coth(double z){ return(cosh(z) / sinh(z)); }

//==============================================================================
/// Returns secant of angle in radians.
//==============================================================================
inline double sec(double z){ return(1.0 / cos(z)); }

//==============================================================================
/// Returns cosecant of input angle in radians.
//==============================================================================
inline double csc(double z){ return(1.0 / sin(z)); }

}

#endif




