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

#ifndef _JSphGpu_
#define _JSphGpu_

#include "Types.h"
#include "JSphTimersGpu.h"
#include "JSph.h"
#include <string>

class JPartsOut;
class JArraysGpu;
class JCellDivGpu;

class JSphGpu : public JSph
{
private:
  JCellDivGpu* CellDiv;

public:
  typedef struct {
    unsigned forcesfluid;
    unsigned forcesbound;
    unsigned forcessps;
    unsigned forcesdem;
    unsigned forcesfluid_zoleft;
    unsigned forcesfluidcorr_zoleft;
    unsigned forcesbound_zoleft;
    unsigned forcesfluid_zoright;
    unsigned forcesfluidcorr_zoright;
    unsigned forcesbound_zoright;
  }StBlockSizes;

protected:
  std::string PtxasFile;      //Fichero con informacion de ptxas.
  StBlockSizes BlockSizes;    //Almacena configuracion de BlockSizes.
  std::string BlockSizesStr;  //Almacena configuracion de BlockSizes en texto.

  //-Vars. con informacion del hardware GPU.
  int GpuSelect;          ///<Gpu seleccionada (-1:sin seleccion).
  std::string GpuName;    ///<Nombre de la gpu seleccionada.
  size_t GpuGlobalMem;    ///<Tamaño de memoria global en bytes.
  unsigned GpuSharedMem;  ///<Tamaño de memoria shared por bloque en bytes.
  unsigned GpuCompute;    ///<Compute capability: 10,11,12,20

  std::string RunMode;     //-Almacena modo de ejecucion (simetria,openmp,balanceo,...).

  //-Numero de particulas del dominio.
  unsigned Np;     //-Numero total de particulas (incluidas las duplicadas periodicas).
  unsigned Npb;    //-Numero de particulas contorno (incluidas las contorno periodicas).
  unsigned NpbOk;  //-Numero de particulas contorno cerca del fluido (incluidas las contorno periodicas).

  unsigned NpfPer;   //-Numero de particulas fluidas-floating periodicas.
  unsigned NpbPer;   //-Numero de particulas contorno periodicas.
  unsigned NpfPerM1; //-Numero de particulas fluidas-floating periodicas (valores anteriores).
  unsigned NpbPerM1; //-Numero de particulas contorno periodicas (valores anteriores).

  bool WithFloating;
  bool BoundChanged; //-Indica si el contorno seleccionado a cambiado desde el ultimo divide.

  unsigned CpuParticlesSize; //-Numero de particulas para las cuales se reservo memoria en cpu.
  llong MemCpuParticles;     //-Mermoria reservada para vectores de datos de particulas.

  //-Variables con datos de las particulas para ejecucion (size=ParticlesSize).
  unsigned *Idp;   //-Identificador de particula.
  word *Code;      //-Indica el grupo de las particulas y otras marcas especiales.
  unsigned *Dcell; //-Celda dentro de DomCells codificada con DomCellCode.
  tdouble2 *Posxy;
  double *Posz;
  tfloat4 *Velrhop;

  //-Variables auxiliares para conversion (size=ParticlesSize).
  tdouble3 *AuxPos;
  tfloat3 *AuxVel; 
  float *AuxRhop;

  unsigned GpuParticlesSize;  //-Numero de particulas para las cuales se reservo memoria en gpu.
  llong MemGpuParticles;      //-Mermoria reservada para vectores de datos de particulas.
  llong MemGpuFixed;          //-Mermoria reservada en AllocGpuMemoryFixed.

  //-Posicion de particula segun id para motion.
  unsigned *RidpMoveg;//-Solo para boundary moving particles [CaseNmoving] y cuando CaseNmoving!=0 

  //-Lista de arrays en Gpu para particulas.
  JArraysGpu* ArraysGpu;
  //-Variables con datos de las particulas para ejecucion (size=ParticlesSize).
  unsigned *Idpg;   //-Identificador de particula.
  word *Codeg;      //-Indica el grupo de las particulas y otras marcas especiales.
  unsigned *Dcellg; //-Celda dentro de DomCells codificada con DomCellCode.
  double2 *Posxyg;
  double *Poszg;
  float4 *Velrhopg;
    
  //-Vars. para compute step: VERLET
  float4 *VelrhopM1g;  //-Verlet: para guardar valores anteriores
  int VerletStep;

  //-Vars. para compute step: SYMPLECTIC
  double2 *PosxyPreg;  //-Sympletic: para guardar valores en predictor
  double *PoszPreg;
  float4 *VelrhopPreg;
  double DtPre;   

  //-Variables for floating bodies.
  unsigned *FtRidpg;         ///<Identifier to access to the particles of the floating object [CaseNfloat] in GPU.
  float *FtoMasspg;          ///<Mass of the particle for each floating body [FtCount] in GPU (used in interaction forces).

  float4 *FtoDatag;    //-Datos constantes de floatings {pini_u,np_u,radius_f,mass_f} [FtCount]  //__device__ int __float_as_int    (   float   x    )  //__device__ float __int_as_float   (   int     x    )  
  float3 *FtoForcesg;  //-Almacena fuerzas de floatings {face_f3,fomegavel_f3} equivalente a JSphCpu::FtoForces [FtCount].

  double3 *FtoCenterg; //-Mantiene centro de floating. [FtCount] 
  float3 *FtoVelg;     //-Mantiene vel de floating. [FtCount] 
  float3 *FtoOmegag;   //-Mantiene omega de floating. [FtCount] 

  StFtoForces *FtoForces; //-Almacena fuerzas de floatings en CPU [FtCount].
  tdouble3 *FtoCenter;    //-Almacena centro de floating en CPU. [FtCount] 

  //-Variables for DEM. (DEM)
  float4 *DemDatag;          ///<Data of the object {mass, (1-poisson^2)/young, kfric, restitu} in GPU [DemObjsSize].

  //-Vars. para computo de fuerzas.
  float4 *PsPospressg; //-Posicion y press para interaccion Pos-Simple. press=cteb*(powf(rhop*OVERRHOPCERO,gamma)-1.0f);

  float *ViscDtg;
  float3 *Aceg;      //-Acumula fuerzas de interaccion
  float *Arg; 
  float *Deltag;     //-Acumula ajuste de Delta-SPH con DELTA_DynamicExt

  float3 *ShiftPosg;    //-Particle displacement using Shifting.
  float *ShiftDetectg;  //-Used to detect free surface with Shifting.

  double VelMax;      ///<Maximum value of Vel[] sqrt(vel.x^2 + vel.y^2 + vel.z^2) computed in PreInteraction_Forces().
  double AceMax;      ///<Maximum value of Ace[] (ace.x^2 + ace.y^2 + ace.z^2) computed in Interaction_Forces().
  float ViscDtMax;    ///<Valor maximo de ViscDt calculado en Interaction_Forces().

  //-Variables for Laminar+SPS viscosity.  
  tsymatrix3f *SpsTaug;       ///<SPS sub-particle stress tensor.
  tsymatrix3f *SpsGradvelg;   ///<Velocity gradients.


  TimersGpu Timers;


  void InitVars();
  void RunExceptionCuda(const std::string &method,const std::string &msg,cudaError_t error);
  void CheckCudaError(const std::string &method,const std::string &msg);

  void FreeGpuMemoryFixed();
  void AllocGpuMemoryFixed();
  void FreeCpuMemoryParticles();
  void AllocCpuMemoryParticles(unsigned np);
  void FreeGpuMemoryParticles();
  void AllocGpuMemoryParticles(unsigned np,float over);

  void ResizeGpuMemoryParticles(unsigned np);
  void ReserveBasicArraysGpu();

  template<class T> T* TSaveArrayGpu(unsigned np,const T *datasrc)const;
  word*        SaveArrayGpu(unsigned np,const word        *datasrc)const{ return(TSaveArrayGpu<word>       (np,datasrc)); }
  unsigned*    SaveArrayGpu(unsigned np,const unsigned    *datasrc)const{ return(TSaveArrayGpu<unsigned>   (np,datasrc)); }
  float*       SaveArrayGpu(unsigned np,const float       *datasrc)const{ return(TSaveArrayGpu<float>      (np,datasrc)); }
  float4*      SaveArrayGpu(unsigned np,const float4      *datasrc)const{ return(TSaveArrayGpu<float4>     (np,datasrc)); }
  double*      SaveArrayGpu(unsigned np,const double      *datasrc)const{ return(TSaveArrayGpu<double>     (np,datasrc)); }
  double2*     SaveArrayGpu(unsigned np,const double2     *datasrc)const{ return(TSaveArrayGpu<double2>    (np,datasrc)); }
  tsymatrix3f* SaveArrayGpu(unsigned np,const tsymatrix3f *datasrc)const{ return(TSaveArrayGpu<tsymatrix3f>(np,datasrc)); }
  unsigned*    SaveArrayGpu_Uint(unsigned np,const unsigned *datasrc)const;
  template<class T> void TRestoreArrayGpu(unsigned np,T *data,T *datanew)const;
  void RestoreArrayGpu(unsigned np,word        *data,word        *datanew)const{ TRestoreArrayGpu<word>       (np,data,datanew); }
  void RestoreArrayGpu(unsigned np,unsigned    *data,unsigned    *datanew)const{ TRestoreArrayGpu<unsigned>   (np,data,datanew); }
  void RestoreArrayGpu(unsigned np,float       *data,float       *datanew)const{ TRestoreArrayGpu<float>      (np,data,datanew); }
  void RestoreArrayGpu(unsigned np,float4      *data,float4      *datanew)const{ TRestoreArrayGpu<float4>     (np,data,datanew); }
  void RestoreArrayGpu(unsigned np,double      *data,double      *datanew)const{ TRestoreArrayGpu<double>     (np,data,datanew); }
  void RestoreArrayGpu(unsigned np,double2     *data,double2     *datanew)const{ TRestoreArrayGpu<double2>    (np,data,datanew); }
  void RestoreArrayGpu(unsigned np,tsymatrix3f *data,tsymatrix3f *datanew)const{ TRestoreArrayGpu<tsymatrix3f>(np,data,datanew); }
  void RestoreArrayGpu_Uint(unsigned np,unsigned *data,unsigned *datanew)const;

  llong GetAllocMemoryCpu()const;
  llong GetAllocMemoryGpu()const;
  void PrintAllocMemory(llong mcpu,llong mgpu)const;

  void ConstantDataUp();
  void ParticlesDataUp(unsigned n);
  unsigned ParticlesDataDown(unsigned n,unsigned pini,bool code,bool cellorderdecode,bool onlynormal);
  
  void SelecDevice(int gpuid);
  static unsigned OptimizeBlockSize(unsigned compute,unsigned nreg);
  unsigned BlockSizeConfig(const std::string& opname,unsigned compute,tuint2 data);
  void ConfigBlockSizes(bool usezone,bool useperi);

  void ConfigRunMode(std::string preinfo);
  void ConfigCellDiv(JCellDivGpu* celldiv){ CellDiv=celldiv; }
  void InitFloating();
  void InitRun();

  void PreInteractionVars_Forces(TpInter tinter,unsigned np,unsigned npb);
  void PreInteraction_Forces(TpInter tinter);
  void PosInteraction_Forces();
  
  void ComputeVerlet(double dt);
  void ComputeSymplecticPre(double dt);
  void ComputeSymplecticCorr(double dt);
  double DtVariable(bool final);
  void RunShifting(double dt);

  void RunMotion(double stepdt);

  void ShowTimers(bool onlyfile=false);
  void GetTimersInfo(std::string &hinfo,std::string &dinfo)const;
  unsigned TimerGetCount()const{ return(TmgGetCount()); }
  bool TimerIsActive(unsigned ct)const{ return(TmgIsActive(Timers,(CsTypeTimerGPU)ct)); }
  float TimerGetValue(unsigned ct)const{ return(TmgGetValue(Timers,(CsTypeTimerGPU)ct)); }
  const double* TimerGetPtrValue(unsigned ct)const{ return(TmgGetPtrValue(Timers,(CsTypeTimerGPU)ct)); }
  std::string TimerGetName(unsigned ct)const{ return(TmgGetName((CsTypeTimerGPU)ct)); }
  std::string TimerToText(unsigned ct)const{ return(JSph::TimerToText(TimerGetName(ct),TimerGetValue(ct))); }

public:
  JSphGpu(bool withmpi);
  ~JSphGpu();

};

#endif


