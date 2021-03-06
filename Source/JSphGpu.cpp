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

#include "JSphGpu.h"
#include "JSphGpu_ker.h"
#include "JPtxasInfo.h"
#include "JCellDivGpu.h"
#include "JPartFloatBi4.h"
#include "Functions.h"
#include "JSphMotion.h"
#include "JArraysGpu.h"
#include "JSphDtFixed.h"
#include "JSaveDt.h"
#include "JWaveGen.h"
#include "JXml.h"

#include "JBuffer.h"
#include "JFormatFiles2.h"

using namespace std;

//==============================================================================
// Constructor.
//==============================================================================
JSphGpu::JSphGpu(bool withmpi):JSph(false,withmpi){
  ClassName="JSphGpu";
  Idp=NULL; Code=NULL; Dcell=NULL; Posxy=NULL; Posz=NULL; Velrhop=NULL;
  AuxPos=NULL; AuxVel=NULL; AuxRhop=NULL;
  FtoForces=NULL; FtoCenter=NULL;   //-Floatings.
  CellDiv=NULL;
  ArraysGpu=new JArraysGpu;
  InitVars();
  TmgCreation(Timers,false);
}

//==============================================================================
// Destructor.
//==============================================================================
JSphGpu::~JSphGpu(){
  FreeCpuMemoryParticles();
  FreeGpuMemoryParticles();
  FreeGpuMemoryFixed();
  delete ArraysGpu;
  TmgDestruction(Timers);
  cudaDeviceReset();
}

//==============================================================================
// Initialization of variables.
//==============================================================================
void JSphGpu::InitVars(){
  RunMode="";
  memset(&BlockSizes,0,sizeof(StBlockSizes));
  BlockSizesStr="";

  Np=Npb=NpbOk=0;
  NpbPer=NpfPer=0;
  WithFloating=false;

  FreeCpuMemoryParticles();
  Idpg=NULL; Codeg=NULL; Dcellg=NULL; Posxyg=NULL; Poszg=NULL; Velrhopg=NULL;
  VelrhopM1g=NULL;                                 //-Verlet
  PosxyPreg=NULL; PoszPreg=NULL; VelrhopPreg=NULL; //-Symplectic
  PsPospressg=NULL;                                //-Interaccion Pos-Simple.
  SpsTaug=NULL; SpsGradvelg=NULL;                  //-Laminar+SPS. 
  ViscDtg=NULL; 
  Arg=NULL; Aceg=NULL; Deltag=NULL;
  ShiftPosg=NULL; ShiftDetectg=NULL; //-Shifting.
  RidpMoveg=NULL;
  FtRidpg=NULL;    FtoMasspg=NULL;               //-Floatings.
  FtoDatag=NULL;   FtoForcesg=NULL;              //-Calculo de fuerzas en floatings.  
  FtoCenterg=NULL; FtoVelg=NULL; FtoOmegag=NULL; //-Gestion de floatings.
  DemDatag=NULL; //(DEM)
  FreeGpuMemoryParticles();
  FreeGpuMemoryFixed();
}

//==============================================================================
// Lanza excepcion por un error Cuda.
//==============================================================================
void JSphGpu::RunExceptionCuda(const std::string &method,const std::string &msg,cudaError_t error){
  std::string tx=fun::PrintStr("%s (CUDA error: %s).\n",msg.c_str(),cudaGetErrorString(error)); 
  Log->Print(GetExceptionText(method,tx));
  RunException(method,msg);
}

//==============================================================================
// Comprueba error y lanza excepcion si lo hubiera.
//==============================================================================
void JSphGpu::CheckCudaError(const std::string &method,const std::string &msg){
  cudaError_t err=cudaGetLastError();
  if(err!=cudaSuccess)RunExceptionCuda(method,msg,err);
}

//==============================================================================
// Libera memoria fija en Gpu para moving y floating.
//==============================================================================
void JSphGpu::FreeGpuMemoryFixed(){
  MemGpuFixed=0;
  if(RidpMoveg)cudaFree(RidpMoveg);   RidpMoveg=NULL;
  if(FtRidpg)cudaFree(FtRidpg);       FtRidpg=NULL;
  if(FtoMasspg)cudaFree(FtoMasspg);   FtoMasspg=NULL;
  if(FtoDatag)cudaFree(FtoDatag);     FtoDatag=NULL;
  if(FtoForcesg)cudaFree(FtoForcesg); FtoForcesg=NULL;
  if(FtoCenterg)cudaFree(FtoCenterg); FtoCenterg=NULL;
  if(FtoVelg)cudaFree(FtoVelg);       FtoVelg=NULL;
  if(FtoOmegag)cudaFree(FtoOmegag);   FtoOmegag=NULL;
  if(DemDatag)cudaFree(DemDatag);     DemDatag=NULL;
}

//==============================================================================
// Allocates memory for arrays with fixed size (motion and floating bodies).
//==============================================================================
void JSphGpu::AllocGpuMemoryFixed(){
  MemGpuFixed=0;
  //-Allocates memory for moving objects.
  if(CaseNmoving){
    size_t m=sizeof(unsigned)*CaseNmoving;
    cudaMalloc((void**)&RidpMoveg,m);   MemGpuFixed+=m;
  }
  //-Allocates memory for floating bodies.
  if(CaseNfloat){
    size_t m=sizeof(unsigned)*CaseNfloat;
    cudaMalloc((void**)&FtRidpg,m);     MemGpuFixed+=m;
    m=sizeof(float)*FtCount;
    cudaMalloc((void**)&FtoMasspg,m);   MemGpuFixed+=m;
    m=sizeof(float4)*FtCount;
    cudaMalloc((void**)&FtoDatag,m);    MemGpuFixed+=m;
    m=sizeof(float3)*2*FtCount;
    cudaMalloc((void**)&FtoForcesg,m);  MemGpuFixed+=m;
    m=sizeof(double3)*FtCount;
    cudaMalloc((void**)&FtoCenterg,m);  MemGpuFixed+=m;
    m=sizeof(float3)*CaseNfloat;
    cudaMalloc((void**)&FtoVelg,m); MemGpuFixed+=m;
    m=sizeof(float3)*FtCount;
    cudaMalloc((void**)&FtoOmegag,m); MemGpuFixed+=m;
  }
  if(UseDEM){ //(DEM)
    size_t m=sizeof(float4)*DemObjsSize;
    cudaMalloc((void**)&DemDatag,m);    MemGpuFixed+=m;
  }
}

//==============================================================================
// Libera memoria para datos principales de particulas.
//==============================================================================
void JSphGpu::FreeCpuMemoryParticles(){
  CpuParticlesSize=0;
  MemCpuParticles=0;
  delete[] Idp;        Idp=NULL;
  delete[] Code;       Code=NULL;
  delete[] Dcell;      Dcell=NULL;
  delete[] Posxy;      Posxy=NULL;
  delete[] Posz;       Posz=NULL;
  delete[] Velrhop;    Velrhop=NULL;
  delete[] AuxPos;     AuxPos=NULL;
  delete[] AuxVel;     AuxVel=NULL;
  delete[] AuxRhop;    AuxRhop=NULL;
  delete[] FtoForces;  FtoForces=NULL;
  delete[] FtoCenter;  FtoCenter=NULL;
}

//==============================================================================
// Reserva memoria para datos principales de particulas.
//==============================================================================
void JSphGpu::AllocCpuMemoryParticles(unsigned np){
  const char* met="AllocCpuMemoryParticles";
  FreeCpuMemoryParticles();
  CpuParticlesSize=np;
  if(np>0){
    try{
      Idp=new unsigned[np];      MemCpuParticles+=sizeof(unsigned)*np;
      Code=new word[np];         MemCpuParticles+=sizeof(word)*np;
      Dcell=new unsigned[np];    MemCpuParticles+=sizeof(unsigned)*np;
      Posxy=new tdouble2[np];    MemCpuParticles+=sizeof(tdouble2)*np;
      Posz=new double[np];       MemCpuParticles+=sizeof(double)*np;
      Velrhop=new tfloat4[np];   MemCpuParticles+=sizeof(tfloat4)*np;
      AuxPos=new tdouble3[np];   MemCpuParticles+=sizeof(tdouble3)*np; 
      AuxVel=new tfloat3[np];    MemCpuParticles+=sizeof(tfloat3)*np;
      AuxRhop=new float[np];     MemCpuParticles+=sizeof(float)*np;
      //-Memoria auxiliar para floatings.
      FtoForces=new StFtoForces[FtCount];  MemCpuParticles+=sizeof(StFtoForces)*FtCount;
      FtoCenter=new tdouble3[FtCount];     MemCpuParticles+=sizeof(tdouble3)*FtCount;
    }
    catch(const std::bad_alloc){
      RunException(met,fun::PrintStr("Could not allocate the requested memory (np=%u).",np));
    }
  }
}

//==============================================================================
// Libera memoria en Gpu para particulas.
//==============================================================================
void JSphGpu::FreeGpuMemoryParticles(){
  GpuParticlesSize=0;
  MemGpuParticles=0;
  ArraysGpu->Reset();
}

//==============================================================================
// Reserva memoria en Gpu para las particulas. 
//==============================================================================
void JSphGpu::AllocGpuMemoryParticles(unsigned np,float over){
  const char* met="AllocGpuMemoryParticles";
  FreeGpuMemoryParticles();
  //-Calcula numero de particulas para las que se reserva memoria.
  const unsigned np2=(over>0? unsigned(over*np): np);
  GpuParticlesSize=np2;
  //-Calcula cuantos arrays.
  ArraysGpu->SetArraySize(np2);
  ArraysGpu->AddArrayCount(JArraysGpu::SIZE_2B,2);  //-code,code2
  ArraysGpu->AddArrayCount(JArraysGpu::SIZE_4B,4);  //-idp,ar,viscdt,dcell
  if(TDeltaSph==DELTA_DynamicExt)ArraysGpu->AddArrayCount(JArraysGpu::SIZE_4B,1);  //-delta
  ArraysGpu->AddArrayCount(JArraysGpu::SIZE_12B,1); //-ace
  ArraysGpu->AddArrayCount(JArraysGpu::SIZE_16B,4); //-velrhop,posxy
  ArraysGpu->AddArrayCount(JArraysGpu::SIZE_8B,2);  //-posz
  if(TStep==STEP_Verlet){
    ArraysGpu->AddArrayCount(JArraysGpu::SIZE_16B,1); //-velrhopm1
  }
  else if(TStep==STEP_Symplectic){
    ArraysGpu->AddArrayCount(JArraysGpu::SIZE_8B,1);  //-poszpre
    ArraysGpu->AddArrayCount(JArraysGpu::SIZE_16B,2); //-posxypre,velrhoppre
  }
  if(TVisco==VISCO_LaminarSPS){     
    ArraysGpu->AddArrayCount(JArraysGpu::SIZE_24B,2); //-SpsTau,SpsGradvel
  }
  if(CaseNfloat){
    ArraysGpu->AddArrayCount(JArraysGpu::SIZE_4B,4);  //-FtMasspg
  }
  if(TShifting!=SHIFT_None){
    ArraysGpu->AddArrayCount(JArraysGpu::SIZE_12B,1); //-shiftpos
    if(ShiftTFS)ArraysGpu->AddArrayCount(JArraysGpu::SIZE_4B,1); //-shiftdetectc
  }
  if(RenCorrection){
    ArraysGpu->AddArrayCount(JArraysGpu::SIZE_4B,1); //-presskf
  }
  //-Muestra la memoria reservada.
  MemGpuParticles=ArraysGpu->GetAllocMemoryGpu();
  PrintSizeNp(np2,MemGpuParticles);
  CheckCudaError(met,"Failed GPU memory allocation.");
}

//==============================================================================
// Resizes space in GPU memory for particles.
//==============================================================================
void JSphGpu::ResizeGpuMemoryParticles(unsigned npnew){
  //-Saves current data from GPU.
  unsigned    *idp       =SaveArrayGpu(Np,Idpg);
  word        *code      =SaveArrayGpu(Np,Codeg);
  unsigned    *dcell     =SaveArrayGpu(Np,Dcellg);
  double2     *posxy     =SaveArrayGpu(Np,Posxyg);
  double      *posz      =SaveArrayGpu(Np,Poszg);
  float4      *velrhop   =SaveArrayGpu(Np,Velrhopg);
  float4      *velrhopm1 =SaveArrayGpu(Np,VelrhopM1g);
  double2     *posxypre  =SaveArrayGpu(Np,PosxyPreg);
  double      *poszpre   =SaveArrayGpu(Np,PoszPreg);
  float4      *velrhoppre=SaveArrayGpu(Np,VelrhopPreg);
  tsymatrix3f *spstau    =SaveArrayGpu(Np,SpsTaug);
  //-Frees pointers.
  ArraysGpu->Free(Idpg);
  ArraysGpu->Free(Codeg);
  ArraysGpu->Free(Dcellg);
  ArraysGpu->Free(Posxyg);
  ArraysGpu->Free(Poszg);
  ArraysGpu->Free(Velrhopg);
  ArraysGpu->Free(VelrhopM1g);
  ArraysGpu->Free(PosxyPreg);
  ArraysGpu->Free(PoszPreg);
  ArraysGpu->Free(VelrhopPreg);
  ArraysGpu->Free(SpsTaug);
  //-Resizes GPU memory allocation.
  //sprintf(Cad,"++>ResizeGpuMemoryParticles> GpuParticlesSize:%u npnew:%u",GpuParticlesSize,npnew); Log->PrintDbg(Cad);
  const double mbparticle=(double(MemGpuParticles)/(1024*1024))/GpuParticlesSize; //-MB por particula.
  Log->Printf("**JSphGpu: Requesting gpu memory for %u particles: %.1f MB.",npnew,mbparticle*npnew);
  ArraysGpu->SetArraySize(npnew);
  //-Reserve pointers.
  Idpg    =ArraysGpu->ReserveUint();
  Codeg   =ArraysGpu->ReserveWord();
  Dcellg  =ArraysGpu->ReserveUint();
  Posxyg  =ArraysGpu->ReserveDouble2();
  Poszg   =ArraysGpu->ReserveDouble();
  Velrhopg=ArraysGpu->ReserveFloat4();
  if(velrhopm1) VelrhopM1g =ArraysGpu->ReserveFloat4();
  if(posxypre)  PosxyPreg  =ArraysGpu->ReserveDouble2();
  if(poszpre)   PoszPreg   =ArraysGpu->ReserveDouble();
  if(velrhoppre)VelrhopPreg=ArraysGpu->ReserveFloat4();
  if(spstau)    SpsTaug    =ArraysGpu->ReserveSymatrix3f();
  //-Restore data in GPU memory.
  RestoreArrayGpu(Np,idp,Idpg);
  RestoreArrayGpu(Np,code,Codeg);
  RestoreArrayGpu(Np,dcell,Dcellg);
  RestoreArrayGpu(Np,posxy,Posxyg);
  RestoreArrayGpu(Np,posz,Poszg);
  RestoreArrayGpu(Np,velrhop,Velrhopg);
  RestoreArrayGpu(Np,velrhopm1,VelrhopM1g);
  RestoreArrayGpu(Np,posxypre,PosxyPreg);
  RestoreArrayGpu(Np,poszpre,PoszPreg);
  RestoreArrayGpu(Np,velrhoppre,VelrhopPreg);
  RestoreArrayGpu(Np,spstau,SpsTaug);
  //-Updates values.
  GpuParticlesSize=npnew;
  MemGpuParticles=ArraysGpu->GetAllocMemoryGpu();
}

//==============================================================================
// Saves a GPU array in CPU memory. 
//==============================================================================
template<class T> T* JSphGpu::TSaveArrayGpu(unsigned np,const T *datasrc)const{
  T *data=NULL;
  if(datasrc){
    try{
      data=new T[np];
    }
    catch(const std::bad_alloc){
      RunException("TSaveArrayGpu","Could not allocate the requested memory.");
    }
    cudaMemcpy(data,datasrc,sizeof(T)*np,cudaMemcpyDeviceToHost);
  }
  return(data);
}
//==============================================================================
unsigned* JSphGpu::SaveArrayGpu_Uint(unsigned np,const unsigned *datasrc)const{
  unsigned *data=NULL;
  if(datasrc){
    try{
      data=new unsigned[np];
    }
    catch(const std::bad_alloc){
      RunException("SaveArrayGpu_Uint","Could not allocate the requested memory.");
    }
    cudaMemcpy(data,datasrc,sizeof(unsigned)*np,cudaMemcpyDeviceToHost);
  }
  return(data);
}

//==============================================================================
// Restores a GPU array from CPU memory. 
//==============================================================================
template<class T> void JSphGpu::TRestoreArrayGpu(unsigned np,T *data,T *datanew)const{
  if(data&&datanew)cudaMemcpy(datanew,data,sizeof(T)*np,cudaMemcpyHostToDevice);
  delete[] data;
}
//==============================================================================
void JSphGpu::RestoreArrayGpu_Uint(unsigned np,unsigned *data,unsigned *datanew)const{
  if(data&&datanew)cudaMemcpy(datanew,data,sizeof(unsigned)*np,cudaMemcpyHostToDevice);
  delete[] data;
}

//==============================================================================
// Arrays para datos basicos de las particulas. 
//==============================================================================
void JSphGpu::ReserveBasicArraysGpu(){
  Idpg=ArraysGpu->ReserveUint();
  Codeg=ArraysGpu->ReserveWord();
  Dcellg=ArraysGpu->ReserveUint();
  Posxyg=ArraysGpu->ReserveDouble2();
  Poszg=ArraysGpu->ReserveDouble();
  Velrhopg=ArraysGpu->ReserveFloat4();
  if(TStep==STEP_Verlet)VelrhopM1g=ArraysGpu->ReserveFloat4();
  if(TVisco==VISCO_LaminarSPS)SpsTaug=ArraysGpu->ReserveSymatrix3f();
}

//==============================================================================
// Devuelve la memoria reservada en cpu.
//==============================================================================
llong JSphGpu::GetAllocMemoryCpu()const{  
  llong s=JSph::GetAllocMemoryCpu();
  //Reservada en AllocMemoryParticles()
  s+=MemCpuParticles;
  //Reservada en otros objetos
  return(s);
}

//==============================================================================
// Devuelve la memoria reservada en gpu.
//==============================================================================
llong JSphGpu::GetAllocMemoryGpu()const{  
  long long s=0;
  //Reservada en AllocGpuMemoryParticles()
  s+=MemGpuParticles;
  //Reservada en AllocGpuMemoryFixed()
  s+=MemGpuFixed;
  //Reservada en otros objetos
  return(s);
}

//==============================================================================
// Visualiza la memoria reservada
//==============================================================================
void JSphGpu::PrintAllocMemory(llong mcpu,llong mgpu)const{
  Log->Printf("Allocated memory in CPU: %lld (%.2f MB)",mcpu,double(mcpu)/(1024*1024));
  Log->Printf("Allocated memory in GPU: %lld (%.2f MB)",mgpu,double(mgpu)/(1024*1024));
}

//==============================================================================
// Copia constantes a la GPU.
//==============================================================================
void JSphGpu::ConstantDataUp(){
  StCteInteraction ctes;
  memset(&ctes,0,sizeof(StCteInteraction));
  ctes.nbound=CaseNbound;
  ctes.massb=MassBound; ctes.massf=MassFluid;
  ctes.fourh2=Fourh2; ctes.h=H;
  ctes.awen=Awen; ctes.bwen=Bwen;
  ctes.cs0=float(Cs0); ctes.eta2=Eta2;
  ctes.delta2h=Delta2H;
  ctes.scell=Scell; ctes.dosh=Dosh; ctes.dp=float(Dp);
  ctes.cteb=CteB; ctes.gamma=Gamma;
  ctes.rhopzero=RhopZero;
  ctes.ovrhopzero=1.f/RhopZero;
  ctes.movlimit=MovLimit;
  ctes.maprealposminx=MapRealPosMin.x; ctes.maprealposminy=MapRealPosMin.y; ctes.maprealposminz=MapRealPosMin.z;
  ctes.maprealsizex=MapRealSize.x; ctes.maprealsizey=MapRealSize.y; ctes.maprealsizez=MapRealSize.z;
  ctes.periactive=PeriActive;
  ctes.xperincx=PeriXinc.x; ctes.xperincy=PeriXinc.y; ctes.xperincz=PeriXinc.z;
  ctes.yperincx=PeriYinc.x; ctes.yperincy=PeriYinc.y; ctes.yperincz=PeriYinc.z;
  ctes.zperincx=PeriZinc.x; ctes.zperincy=PeriZinc.y; ctes.zperincz=PeriZinc.z;
  ctes.cellcode=DomCellCode;
  ctes.domposminx=DomPosMin.x; ctes.domposminy=DomPosMin.y; ctes.domposminz=DomPosMin.z;
  cusph::CteInteractionUp(&ctes);
  CheckCudaError("ConstantDataUp","Failed copying constants to GPU.");
}

//==============================================================================
// Sube datos de particulas a la GPU.
//==============================================================================
void JSphGpu::ParticlesDataUp(unsigned n){
  cudaMemcpy(Idpg,Idp,sizeof(unsigned)*n,cudaMemcpyHostToDevice);
  cudaMemcpy(Codeg,Code,sizeof(word)*n,cudaMemcpyHostToDevice);
  cudaMemcpy(Dcellg,Dcell,sizeof(unsigned)*n,cudaMemcpyHostToDevice);
  cudaMemcpy(Posxyg,Posxy,sizeof(double2)*n,cudaMemcpyHostToDevice);
  cudaMemcpy(Poszg,Posz,sizeof(double)*n,cudaMemcpyHostToDevice);
  cudaMemcpy(Velrhopg,Velrhop,sizeof(float4)*n,cudaMemcpyHostToDevice);
  CheckCudaError("ParticlesDataUp","Failed copying data to GPU.");
}

//==============================================================================
// Recupera datos de particulas de la GPU y devuelve el numero de particulas que
// sera menor que n si se eliminaron las periodicas.
// - code: Recupera datos de Codeg.
// - cellorderdecode: Reordena componentes de pos y vel segun CellOrder.
// - onlynormal: Solo se queda con las normales, elimina las particulas periodicas.
//==============================================================================
unsigned JSphGpu::ParticlesDataDown(unsigned n,unsigned pini,bool code,bool cellorderdecode,bool onlynormal){
  unsigned num=n;
  cudaMemcpy(Idp,Idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  cudaMemcpy(Posxy,Posxyg+pini,sizeof(double2)*n,cudaMemcpyDeviceToHost);
  cudaMemcpy(Posz,Poszg+pini,sizeof(double)*n,cudaMemcpyDeviceToHost);
  cudaMemcpy(Velrhop,Velrhopg+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
  if(code || onlynormal)cudaMemcpy(Code,Codeg+pini,sizeof(word)*n,cudaMemcpyDeviceToHost);
  CheckCudaError("ParticlesDataDown","Failed copying data from GPU.");
  //-Elimina particulas no normales (periodicas y otras).
  if(onlynormal){
    unsigned ndel=0;
    for(unsigned p=0;p<n;p++){
      bool normal=(CODE_GetSpecialValue(Code[p])==CODE_NORMAL);
      if(ndel && normal){
        Idp[p-ndel]    =Idp[p];
        Posxy[p-ndel]  =Posxy[p];
        Posz[p-ndel]   =Posz[p];
        Velrhop[p-ndel]=Velrhop[p];
        Code[p-ndel]   =Code[p];
      }
      if(!normal)ndel++;
    }
    num-=ndel;
  }
  //-Convierte datos a formato simple.
  for(unsigned p=0;p<n;p++){
    AuxPos[p]=TDouble3(Posxy[p].x,Posxy[p].y,Posz[p]);
    AuxVel[p]=TFloat3(Velrhop[p].x,Velrhop[p].y,Velrhop[p].z);
    AuxRhop[p]=Velrhop[p].w;
  }
  if(cellorderdecode)DecodeCellOrder(n,AuxPos,AuxVel);
  return(num);
}

//==============================================================================
// Inicializa dispositivo cuda 
//==============================================================================
void JSphGpu::SelecDevice(int gpuid){
  const char* met="SelecDevice";
  Log->Print("[Select CUDA Device]");
  GpuSelect=-1;
  int devcount;
  cudaGetDeviceCount(&devcount);
  CheckCudaError(met,"Failed getting devices info.");
  for(int dev=0;dev<devcount;dev++){
    cudaDeviceProp devp;
    cudaGetDeviceProperties(&devp,dev);
    Log->Printf("Device %d: \"%s\"",dev,devp.name);
    Log->Printf("  Compute capability:        %d.%d",devp.major,devp.minor);
    int corebymp=(devp.major==1? 8: (devp.major==2? (devp.minor==0? 32: 48): (devp.major==3? 192: -1)));
    Log->Printf("  Multiprocessors:           %d (%d cores)",devp.multiProcessorCount,devp.multiProcessorCount*corebymp);
    Log->Printf("  Memory global:             %d MB",int(devp.totalGlobalMem/(1024*1024)));
    Log->Printf("  Clock rate:                %.2f GHz",devp.clockRate*1e-6f);
    #if CUDART_VERSION >= 2020
    Log->Printf("  Run time limit on kernels: %s",(devp.kernelExecTimeoutEnabled? "Yes": "No"));
    #endif
    #if CUDART_VERSION >= 3010
    Log->Printf("  ECC support enabled:       %s",(devp.ECCEnabled? "Yes": "No"));
    #endif
  }
  Log->Print("");
  if(devcount){
    if(gpuid>=0)cudaSetDevice(gpuid);
    else{
      unsigned *ptr=NULL;
      cudaMalloc((void**)&ptr,sizeof(unsigned)*100);
      cudaFree(ptr);
    }
    cudaDeviceProp devp;
    int dev;
    cudaGetDevice(&dev);
    cudaGetDeviceProperties(&devp,dev);
    GpuSelect=dev;
    GpuName=devp.name;
    GpuGlobalMem=devp.totalGlobalMem;
    GpuSharedMem=int(devp.sharedMemPerBlock);
    GpuCompute=devp.major*10+devp.minor;
    //-Muestra informacion del hardware seleccionado.
    Log->Print("[GPU Hardware]");
    if(gpuid<0)Hardware=fun::PrintStr("Gpu_%d?=\"%s\"",GpuSelect,GpuName.c_str());
    else Hardware=fun::PrintStr("Gpu_%d=\"%s\"",GpuSelect,GpuName.c_str());

    if(gpuid<0)Log->Printf("Device default: %d  \"%s\"",GpuSelect,GpuName.c_str());
    else Log->Printf("Device selected: %d  \"%s\"",GpuSelect,GpuName.c_str());
    Log->Printf("Compute capability: %.1f",float(GpuCompute)/10);
    Log->Printf("Memory global: %d MB",int(GpuGlobalMem/(1024*1024)));
    Log->Printf("Memory shared: %u Bytes",GpuSharedMem);
  }
  else RunException(met,"No hay dispositivos CUDA disponibles.");
}

//==============================================================================
// Devuelve el tama�o optimo de bloque segun registros del kernel y compute.
//==============================================================================
unsigned JSphGpu::OptimizeBlockSize(unsigned compute,unsigned nreg){
  //printf("compute:%d\n",compute);
  if(compute>=30){
    if(nreg<=32)return(256);       // 1-32 -> 128:100%  256:100%  512:100%
    else if(nreg<=40)return(256);  //33-40 -> 128:75%  256:75%  384:75%  512:75%
    else if(nreg<=48)return(256);  //41-48 -> 128:63%  256:63%
    else if(nreg<=56)return(128);  //49-56 -> 128:56%  256:50%  384:56%
    else if(nreg<=63)return(256);  //49-63 -> 128:50%  256:50%  512:50%
    //-For Compute capability 3.5
    else if(nreg<=64)return(256);  //64      -> 128:50%  256:50%
    else if(nreg<=72)return(128);  //65-72   -> 128:44%  256:38%
    else if(nreg<=80)return(256);  //73-80   -> 128:38%  256:38%
    else if(nreg<=96)return(128);  //81-96   -> 128:31%  256:25%
    else if(nreg<=128)return(256); //97-128  -> 128:25%  256:25%
    else if(nreg<=168)return(128); //129-168 -> 128:19%  256:13%
    else if(nreg<=255)return(256); //199-255 -> 128:13%  256:13%
    else return(256);              
  }
  else if(compute>=20){
    if(nreg<=20)return(256);       // 1-20 -> 192:100%  256:100%  384:100%
    else if(nreg<=24)return(224);  //21-24 -> 192:88%  224:88%  256:83%  448:88%
    else if(nreg<=26)return(192);  //25-26 -> 192:75%  256:67%  288:75%  416:81%
    else if(nreg<=28)return(192);  //27-28 -> 192:75%  256:67%  288:75%
    else if(nreg<=32)return(256);  //29-32 -> 256:67%
    else if(nreg<=34)return(192);  //33-34 -> 192:63%  256:50%
    else if(nreg<=36)return(128);  //35-36 -> 128:58%  224:58%  256:50%
    else if(nreg<=42)return(128);  //37-42 -> 128:50%  256:50%
    else if(nreg<=44)return(224);  //43-44 -> 192:38%  224:44%  256:33%  352:46%
    else if(nreg<=50)return(128);  //45-50 -> 128:42%  160:42%  256:33%  320:42%
    else if(nreg<=56)return(192);  //51-56 -> 192:38%  256:33%  288:38%
    else if(nreg<=63)return(128);  //57-63 -> 128:33%  256:33%
    else return(128);              
  }
  else if(compute>=12){
    if(nreg<=16)return(256);       // 1-16 -> 128:100%  256:100%
    else if(nreg<=18)return(448);  //17-18 -> 128:75%  192:75%  256:75%  448:88%
    else if(nreg<=20)return(256);  //19-20 -> 192:75%  256:75%  384:75%
    else if(nreg<=21)return(192);  //21    -> 192:75%  256:50%  288:56%  320:63%  352:69%  384:75%
    else if(nreg<=24)return(128);  //22-24 -> 128:63%  192:56%  288:56%  256:50%  320:63%
    else if(nreg<=25)return(320);  //25    -> 192:56%  288:56%  256:50%  320:63%
    else if(nreg<=26)return(192);  //26    -> 192:56%  256:50%
    else if(nreg<=32)return(256);  //27-32 -> 256:50%
    else if(nreg<=36)return(448);  //33-36 -> 192:38%  256:25%  416:41%  448:44%
    else if(nreg<=42)return(192);  //37-42 -> 192:38%  256:25%
    else if(nreg<=51)return(320);  //43-51 -> 256:25%  288:28%  320:31%
    else if(nreg<=64)return(256);  //52-64 -> 128:25%  256:25%
    else return(192);              //65-85 -> 128:13%  192:19%
  }
  else if(compute>=10){
    if(nreg<=10)return(256);       // 1-10 -> 128:100%  192:100%  256:100%  384:100%
    else if(nreg<=12)return(128);  //11-12 -> 128:83%  192:75%  256:67%  320:83%
    else if(nreg<=13)return(192);  //13    -> 128:67%  192:75%  256:67%
    else if(nreg<=16)return(256);  //14-16 -> 128:67%  192:50%  256:67%
    else if(nreg<=18)return(448);  //17-18 -> 128:50%  192:50%  256:33%  384:50%  448:58%
    else if(nreg<=20)return(128);  //19-20 -> 128:50%  192:50%  256:33%  384:50%
    else if(nreg<=21)return(192);  //21    -> 128:33%  192:50%  256:33%  384:50%
    else if(nreg<=24)return(320);  //22-24 -> 64:42%  128:33%  256:33%  320:42%
    else if(nreg<=25)return(320);  //25    -> 64:33%  128:33%  256:33%  320:42%
    else if(nreg<=32)return(256);  //26-32 -> 64:33%  128:33%  256:33%
    else if(nreg<=40)return(192);  //33-40 -> 64:25%  128:17%  192:25%
    else if(nreg<=42)return(192);  //41-42 -> 64:17%  128:17%  192:25%
    else if(nreg<=64)return(128);  //43-64 -> 64:17%  128:17%
    else return(64);               //65-128-> 64:8%
  }
  return(256);
}

//==============================================================================
// Devuelve BlockSize en funcion de registros del kernel.
//==============================================================================
unsigned JSphGpu::BlockSizeConfig(const string& opname,unsigned compute,tuint2 data){
  std::string tx;
  unsigned bsize=256;
  if(data.x){
    bsize=OptimizeBlockSize(compute,data.x);
    if(!data.y)tx=fun::PrintStr("%s=%u (%u regs)",opname.c_str(),bsize,data.x);
    else tx=fun::PrintStr("%s=%u (%u regs + %u bytes)",opname.c_str(),bsize,data.x,data.y);
  }
  else tx=fun::PrintStr("%s=%u (? regs)",opname.c_str(),bsize);
  Log->Print(tx);
  if(!BlockSizesStr.empty())BlockSizesStr=BlockSizesStr+", ";
  BlockSizesStr=BlockSizesStr+tx;
  return(bsize);
}

//==============================================================================
// Configura datos de DeviceContext y DeviceCtes. Devuelve true en caso de error.
//==============================================================================
void JSphGpu::ConfigBlockSizes(bool usezone,bool useperi){
  const char met[]="ConfigBlockSizes";
  //-Obtiene configuracion segun CellMode
  //--------------------------------------
  Log->Print(" ");
  Log->Print(fun::VarStr("PtxasFile",PtxasFile));
  const unsigned smgpu=(GpuCompute<35? (GpuCompute<30? (GpuCompute<20? (GpuCompute<12? 10: 13): 20): 30): 35);
  unsigned smcode=smgpu;//(smgpu==13? 12: smgpu);
  JPtxasInfo pt;
  if(fun::FileExists(PtxasFile)){
    pt.LoadFile(PtxasFile);
    //RunException(met,"Stop...");
    if(smgpu==20&&!pt.CheckSm(20))RunException(met,"Code is not compiled for sm20.");
    if(smgpu==30&&!pt.CheckSm(30)){
      if(!pt.CheckSm(20))RunException(met,"Code is not compiled for sm20 and sm30.");
      else smcode=20;
    }
    if(smgpu==35&&!pt.CheckSm(35)){
      if(!pt.CheckSm(20))RunException(met,"Code is not compiled for sm20 and sm35.");
      else smcode=20;
    }
    Log->Printf("Use code for compute capability %3.1f on hardware %3.1f",float(smcode)/10,float(smgpu)/10);
  }
  else Log->Print("**Without optimization of registers.");
  pt.SaveCsv(DirOut+"ptxas_info.csv");
  BlockSizesStr="";
  if(CellMode==CELLMODE_2H||CellMode==CELLMODE_H){
    const TpFtMode ftmode=(CaseNfloat? (UseDEM? FTMODE_Dem: FTMODE_Sph): FTMODE_None);
    const bool lamsps=(TVisco==VISCO_LaminarSPS);
    const bool shift=(TShifting!=SHIFT_None);
    BlockSizes.forcesbound=BlockSizeConfig("BsForcesBound",smgpu,pt.GetData("cusph_KerInteractionForcesBound",smcode,Psimple,ftmode));
    BlockSizes.forcesfluid=BlockSizeConfig("BsForcesFluid",smgpu,pt.GetData("cusph_KerInteractionForcesFluid",smcode,Psimple,ftmode,lamsps,TDeltaSph,shift));
    if(UseDEM)BlockSizes.forcesdem=BlockSizeConfig("BsForcesDEM",smgpu,pt.GetData("cusph_KerInteractionForcesDem",smcode,Psimple));
  }
  else RunException(met,"CellMode unrecognised.");
  Log->Print(" ");
}

//==============================================================================
// Configura modo de ejecucion en GPU.
//==============================================================================
void JSphGpu::ConfigRunMode(std::string preinfo){
  #ifndef WIN32
    const int len=128; char hname[len];
    gethostname(hname,len);
    if(!preinfo.empty())preinfo=preinfo+", ";
    preinfo=preinfo+"HostName:"+hname;
  #endif
  RunMode=preinfo+RunMode;
  if(Stable)RunMode=string("Stable, ")+RunMode;
  if(Psimple)RunMode=string("Pos-Simple, ")+RunMode;
  else RunMode=string("Pos-Double, ")+RunMode;
  Log->Print(fun::VarStr("RunMode",RunMode));
  if(!RunMode.empty())RunMode=RunMode+", "+BlockSizesStr;
}

//==============================================================================
/// Adjusts variables of particles of floating bodies.
//==============================================================================
void JSphGpu::InitFloating(){
  if(PartBegin){
    JPartFloatBi4Load ftdata;
    ftdata.LoadFile(PartBeginDir);
    //-Comprueba coincidencia de datos constantes.
    for(unsigned cf=0;cf<FtCount;cf++)ftdata.CheckHeadData(cf,FtObjs[cf].mkbound,FtObjs[cf].begin,FtObjs[cf].count,FtObjs[cf].mass);
    //-Carga datos de PART.
    ftdata.LoadPart(PartBegin);
    for(unsigned cf=0;cf<FtCount;cf++){
      FtObjs[cf].center=OrderCodeValue(CellOrder,ftdata.GetPartCenter(cf));
      FtObjs[cf].fvel=OrderCodeValue(CellOrder,ftdata.GetPartFvel(cf));
      FtObjs[cf].fomega=OrderCodeValue(CellOrder,ftdata.GetPartFomega(cf));
      FtObjs[cf].radius=ftdata.GetHeadRadius(cf);
    }
    DemDtForce=ftdata.GetPartDemDtForce();
  }
  //-Copies massp values to GPU.
  {
    float *massp=new float[FtCount];
    for(unsigned cf=0;cf<FtCount;cf++)massp[cf]=FtObjs[cf].massp;
    cudaMemcpy(FtoMasspg,massp,sizeof(float)*FtCount,cudaMemcpyHostToDevice);
    delete[] massp;
  }
  //-Copies floating values to GPU.
  {
    typedef struct{
      unsigned pini;
      unsigned np;
      float radius;
      float mass;
    }stdata;

    stdata *data=new stdata[FtCount];
    tdouble3 *center=new tdouble3[FtCount];
    tfloat3 *vel=new tfloat3[FtCount];
    tfloat3 *omega=new tfloat3[FtCount];
    for(unsigned cf=0;cf<FtCount;cf++){
      data[cf].pini=FtObjs[cf].begin-CaseNpb;
      data[cf].np=FtObjs[cf].count;
      data[cf].radius=FtObjs[cf].radius;
      data[cf].mass=FtObjs[cf].mass;
      //float pini,np;   //-Daba problemas en Linux...
      //*((unsigned*)&pini)=FtObjs[cf].begin-CaseNpb;
      //*((unsigned*)&np  )=FtObjs[cf].count;
      //data[cf]=TFloat4(pini,np,FtObjs[cf].radius,FtObjs[cf].mass);
      center[cf]=FtObjs[cf].center;
      vel[cf]=FtObjs[cf].fvel;
      omega[cf]=FtObjs[cf].fomega;
    }
    cudaMemcpy(FtoDatag  ,data  ,sizeof(float4) *FtCount,cudaMemcpyHostToDevice);
    cudaMemcpy(FtoCenterg,center,sizeof(double3)*FtCount,cudaMemcpyHostToDevice);
    cudaMemcpy(FtoVelg   ,vel   ,sizeof(float3) *FtCount,cudaMemcpyHostToDevice);
    cudaMemcpy(FtoOmegag ,omega ,sizeof(float3) *FtCount,cudaMemcpyHostToDevice);
    delete[] data;
    delete[] center;
    delete[] vel;
    delete[] omega;
  }
  //-Copies data object for DEM to GPU.
  if(UseDEM){ //(DEM)
    float4 *data=new float4[DemObjsSize];
    for(unsigned c=0;c<DemObjsSize;c++){
      data[c].x=DemObjs[c].mass;
      data[c].y=DemObjs[c].tau;
      data[c].z=DemObjs[c].kfric;
      data[c].w=DemObjs[c].restitu;
    }
    cudaMemcpy(DemDatag,data,sizeof(float4)*DemObjsSize,cudaMemcpyHostToDevice);
    delete[] data;
  }
}

//==============================================================================
// Inicializa vectores y variables para la ejecucion.
//==============================================================================
void JSphGpu::InitRun(){
  const char met[]="InitRun";
  WithFloating=(CaseNfloat>0);
  if(TStep==STEP_Verlet){
    cudaMemcpy(VelrhopM1g,Velrhopg,sizeof(float4)*Np,cudaMemcpyDeviceToDevice);
    VerletStep=0;
  }
  else if(TStep==STEP_Symplectic)DtPre=DtIni;
  if(TVisco==VISCO_LaminarSPS)cudaMemset(SpsTaug,0,sizeof(tsymatrix3f)*Np);
  if(UseDEM)DemDtForce=DtIni; //(DEM)
  if(CaseNfloat)InitFloating();

  //-Adjust paramaters to start.
  PartIni=PartBeginFirst;
  TimeStepIni=(!PartIni? 0: PartBeginTimeStep);
  //-Adjust motion for the instant of the loaded PART.
  if(CaseNmoving){
    MotionTimeMod=(!PartIni? PartBeginTimeStep: 0);
    Motion->ProcesTime(0,TimeStepIni+MotionTimeMod);
  }

  //-Uses Inlet information from PART readed.
  if(PartBeginTimeStep && PartBeginTotalNp){
    TotalNp=PartBeginTotalNp;
    IdMax=unsigned(TotalNp-1);
  }

  //-Prepares WaveGen configuration.
  if(WaveGen){
    Log->Printf("\nWave paddles configuration:");
    WaveGen->Init(TimeMax,Gravity,Simulate2D,CellOrder,MassFluid,Dp,Dosh,Scell,Hdiv,DomPosMin,DomRealPosMin,DomRealPosMax);
    WaveGen->VisuConfig(""," ");
  }

  //-Process Special configurations in XML.
  JXml xml; xml.LoadFile(FileXml);
  //-Configuration of SaveDt.
  if(xml.GetNode("case.execution.special.savedt",false)){
    SaveDt=new JSaveDt(Log);
    SaveDt->Config(&xml,"case.execution.special.savedt",TimeMax,TimePart);
    SaveDt->VisuConfig("\nSaveDt configuration:"," ");
  }

  Part=PartIni; Nstep=0; PartNstep=0; PartOut=0;
  TimeStep=TimeStepIni; TimeStepM1=TimeStep;
  if(DtFixed)DtIni=DtFixed->GetDt(TimeStep,DtIni);
  if(TimersStep)TimersStep->SetInitialTime(float(TimeStep));
  CheckCudaError("InitRun","Failed initializing variables for execution.");
}

//==============================================================================
// Prepara variables para interaccion "INTER_Forces" o "INTER_ForcesCorr".
//==============================================================================
void JSphGpu::PreInteractionVars_Forces(TpInter tinter,unsigned ini,unsigned np,unsigned npb){
  //-Inicializa arrays.
  const unsigned npf=np-npb;
  cudaMemset(ViscDtg+ini,0,sizeof(float)*np);                                //ViscDtg[]=0
  cudaMemset(Arg+ini,0,sizeof(float)*np);                                    //Arg[]=0
  if(Deltag)cudaMemset(Deltag+ini,0,sizeof(float)*np);                       //Deltag[]=0
  if(ShiftPosg)cudaMemset(ShiftPosg+ini,0,sizeof(tfloat3)*np);               //ShiftPosg[]=0
  if(ShiftDetectg)cudaMemset(ShiftDetectg+ini,0,sizeof(float)*np);           //ShiftDetectg[]=0
  cudaMemset(Aceg+ini,0,sizeof(tfloat3)*npb);                                //Aceg[]=(0,0,0) para bound
  cusph::InitArray(npf,Aceg+ini+npb,Gravity);                                //Aceg[]=Gravity para fluid
  if(SpsGradvelg)cudaMemset(SpsGradvelg+ini+npb,0,sizeof(tsymatrix3f)*npf);  //SpsGradvelg[]=(0,0,0,0,0,0).
}

//==============================================================================
// Prepara variables para interaccion "INTER_Forces" o "INTER_ForcesCorr".
//==============================================================================
void JSphGpu::PreInteraction_Forces(TpInter tinter){
  TmgStart(Timers,TMG_CfPreForces);
  //-Asigna memoria.
  Arg=ArraysGpu->ReserveFloat();
  Aceg=ArraysGpu->ReserveFloat3();
  ViscDtg=ArraysGpu->ReserveFloat();
  if(TDeltaSph==DELTA_DynamicExt)Deltag=ArraysGpu->ReserveFloat();
  if(TShifting!=SHIFT_None){
    ShiftPosg=ArraysGpu->ReserveFloat3();
    if(ShiftTFS)ShiftDetectg=ArraysGpu->ReserveFloat();
  }   
  if(TVisco==VISCO_LaminarSPS)SpsGradvelg=ArraysGpu->ReserveSymatrix3f();

  //-Prepara datos para interaccion Pos-Simple.
  if(Psimple){
    PsPospressg=ArraysGpu->ReserveFloat4();
    cusph::PreInteractionSimple(Np,Posxyg,Poszg,Velrhopg,PsPospressg,CteB,Gamma);
  }
  //-Inicializa arrays.
  PreInteractionVars_Forces(tinter,0,Np,Npb);

  //-Calcula VelMax: Se incluyen las particulas floatings y no afecta el uso de condiciones periodicas.
  #ifdef DT_ALLPARTICLES
    const unsigned pini=0;
  #else
    const unsigned pini=Npb;
  #endif
  cusph::ComputeVelMod(Np-pini,Velrhopg+pini,ViscDtg);
  float velmax=cusph::ReduMaxFloat(Np-pini,pini,ViscDtg,CellDiv->GetAuxMem(cusph::ReduMaxFloatSize(Np-pini)));
  VelMax=sqrt(velmax);
  cudaMemset(ViscDtg,0,sizeof(float)*Np);           //ViscDtg[]=0
  ViscDtMax=0;
  CheckCudaError("PreInteraction_Forces","Failed calculating VelMax.");
  TmgStop(Timers,TMG_CfPreForces);
}

//==============================================================================
// Libera memoria asignada de ArraysGpu.
//==============================================================================
void JSphGpu::PosInteraction_Forces(){
  //-Libera memoria asignada en PreInteraction_Forces().
  ArraysGpu->Free(Arg);          Arg=NULL;
  ArraysGpu->Free(Aceg);         Aceg=NULL;
  ArraysGpu->Free(ViscDtg);      ViscDtg=NULL;
  ArraysGpu->Free(Deltag);       Deltag=NULL;
  ArraysGpu->Free(ShiftPosg);    ShiftPosg=NULL;
  ArraysGpu->Free(ShiftDetectg); ShiftDetectg=NULL;
  ArraysGpu->Free(PsPospressg);  PsPospressg=NULL;
  ArraysGpu->Free(SpsGradvelg);  SpsGradvelg=NULL;
}

//==============================================================================
// Actualizacion de particulas segun fuerzas y dt usando Verlet.
//==============================================================================
void JSphGpu::ComputeVerlet(double dt){  //pdtedom
  TmgStart(Timers,TMG_SuComputeStep);
  const bool shift=TShifting!=SHIFT_None;
  VerletStep++;
  //-Asigna memoria para calcular el desplazamiento.
  double2 *movxyg=ArraysGpu->ReserveDouble2();
  double *movzg=ArraysGpu->ReserveDouble();
  //-Calcula desplazamiento, velocidad y densidad.
  if(VerletStep<VerletSteps){
    const double twodt=dt+dt;
    cusph::ComputeStepVerlet(WithFloating,shift,Np,Npb,Velrhopg,VelrhopM1g,Arg,Aceg,ShiftPosg,dt,twodt,RhopOutMin,RhopOutMax,Codeg,movxyg,movzg,VelrhopM1g);
  }
  else{
    cusph::ComputeStepVerlet(WithFloating,shift,Np,Npb,Velrhopg,Velrhopg,Arg,Aceg,ShiftPosg,dt,dt,RhopOutMin,RhopOutMax,Codeg,movxyg,movzg,VelrhopM1g);
    VerletStep=0;
  }
  //-Los nuevos valores se calculan en VelrhopM1g.
  swap(Velrhopg,VelrhopM1g);   //-Intercambia Velrhopg y VelrhopM1g.
  //-Aplica desplazamiento a las particulas fluid no periodicas.
  cusph::ComputeStepPos(PeriActive,WithFloating,Np,Npb,movxyg,movzg,Posxyg,Poszg,Dcellg,Codeg);
  //-Libera memoria asignada al desplazamiento.
  ArraysGpu->Free(movxyg);   movxyg=NULL;
  ArraysGpu->Free(movzg);    movzg=NULL;
  TmgStop(Timers,TMG_SuComputeStep);
}

//==============================================================================
// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Predictor.
//==============================================================================
void JSphGpu::ComputeSymplecticPre(double dt){
  TmgStart(Timers,TMG_SuComputeStep);
  const bool shift=TShifting!=SHIFT_None;
  //-Asigna memoria a variables Pre.
  PosxyPreg=ArraysGpu->ReserveDouble2();
  PoszPreg=ArraysGpu->ReserveDouble();
  VelrhopPreg=ArraysGpu->ReserveFloat4();
  //-Cambia datos a variables Pre para calcular nuevos datos.
  swap(PosxyPreg,Posxyg);      //Es decir... PosxyPre[] <= Posxy[]
  swap(PoszPreg,Poszg);        //Es decir... PoszPre[] <= Posz[]
  swap(VelrhopPreg,Velrhopg);  //Es decir... VelrhopPre[] <= Velrhop[]
  //-Asigna memoria para calcular el desplazamiento.
  double2 *movxyg=ArraysGpu->ReserveDouble2();
  double *movzg=ArraysGpu->ReserveDouble();
  //-Calcula desplazamiento, velocidad y densidad.
  const double dt05=dt*.5;
  cusph::ComputeStepSymplecticPre(WithFloating,shift,Np,Npb,VelrhopPreg,Arg,Aceg,ShiftPosg,dt05,RhopOutMin,RhopOutMax,Codeg,movxyg,movzg,Velrhopg);
  //-Aplica desplazamiento a las particulas fluid no periodicas.
  cusph::ComputeStepPos2(PeriActive,WithFloating,Np,Npb,PosxyPreg,PoszPreg,movxyg,movzg,Posxyg,Poszg,Dcellg,Codeg);
  //-Libera memoria asignada al desplazamiento.
  ArraysGpu->Free(movxyg);   movxyg=NULL;
  ArraysGpu->Free(movzg);    movzg=NULL;
  //-Copia posicion anterior del contorno.
  cudaMemcpy(Posxyg,PosxyPreg,sizeof(double2)*Npb,cudaMemcpyDeviceToDevice);
  cudaMemcpy(Poszg,PoszPreg,sizeof(double)*Npb,cudaMemcpyDeviceToDevice);
  TmgStop(Timers,TMG_SuComputeStep);
}

//==============================================================================
// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Corrector.
//==============================================================================
void JSphGpu::ComputeSymplecticCorr(double dt){
  TmgStart(Timers,TMG_SuComputeStep);
  const bool shift=TShifting!=SHIFT_None;
  //-Asigna memoria para calcular el desplazamiento.
  double2 *movxyg=ArraysGpu->ReserveDouble2();
  double *movzg=ArraysGpu->ReserveDouble();
  //-Calcula desplazamiento, velocidad y densidad.
  const double dt05=dt*.5;
  cusph::ComputeStepSymplecticCor(WithFloating,shift,Np,Npb,VelrhopPreg,Arg,Aceg,ShiftPosg,dt05,dt,RhopOutMin,RhopOutMax,Codeg,movxyg,movzg,Velrhopg);
  //-Aplica desplazamiento a las particulas fluid no periodicas.
  cusph::ComputeStepPos2(PeriActive,WithFloating,Np,Npb,PosxyPreg,PoszPreg,movxyg,movzg,Posxyg,Poszg,Dcellg,Codeg);
  //-Libera memoria asignada al desplazamiento.
  ArraysGpu->Free(movxyg);   movxyg=NULL;
  ArraysGpu->Free(movzg);    movzg=NULL;
  //-Libera memoria asignada a variables Pre en ComputeSymplecticPre().
  ArraysGpu->Free(PosxyPreg);    PosxyPreg=NULL;
  ArraysGpu->Free(PoszPreg);     PoszPreg=NULL;
  ArraysGpu->Free(VelrhopPreg);  VelrhopPreg=NULL;
  TmgStop(Timers,TMG_SuComputeStep);
}

//==============================================================================
// Calcula un Dt variable.
//==============================================================================
double JSphGpu::DtVariable(bool final){
  //-dt1 depends on force per unit mass.
  const double acemax=sqrt(double(AceMax));
  const double dt1=(AceMax? (sqrt(double(H)/AceMax)): DBL_MAX); 
  //-dt2 combines the Courant and the viscous time-step controls.
  const double dt2=double(H)/(max(Cs0,VelMax*10.)+double(H)*ViscDtMax);
  //-dt new value of time step.
  double dt=double(CFLnumber)*min(dt1,dt2);
  //Log->Printf("%u>DTT  dt1:%.20f dt2:%.20f dt:%.20f",Nstep,dt1,dt2,dt); 
  if(DtFixed)dt=DtFixed->GetDt(float(TimeStep),float(dt));
  if(dt<double(DtMin)){ dt=double(DtMin); DtModif++; }
  if(SaveDt && final)SaveDt->AddValues(TimeStep,dt,dt1*CFLnumber,dt2*CFLnumber,AceMax,ViscDtMax,VelMax);
  return(dt);
}

//==============================================================================
// Calcula Shifting final para posicion de particulas.
//==============================================================================
void JSphGpu::RunShifting(double dt){
  const double coeftfs=(Simulate2D? 2.0: 3.0)-ShiftTFS;
  cusph::RunShifting(Np,Npb,dt,ShiftCoef,ShiftTFS,coeftfs,Velrhopg,ShiftDetectg,ShiftPosg);
}

//==============================================================================
// Procesa movimiento de boundary particles
//==============================================================================
void JSphGpu::RunMotion(double stepdt){
  const char met[]="RunMotion";
  TmgStart(Timers,TMG_SuMotion);

  unsigned nmove=0;
  if(Motion->ProcesTime(TimeStep+MotionTimeMod,stepdt)){
    nmove=Motion->GetMovCount();
    //{ char cad[256]; sprintf(cad,"----RunMotion[%u]>  nmove:%u",Nstep,nmove); Log->Print(cad); }
    if(nmove){
      cusph::CalcRidp(PeriActive!=0,Npb,0,CaseNfixed,CaseNfixed+CaseNmoving,Codeg,Idpg,RidpMoveg);
      //-Movimiento de particulas boundary
      for(unsigned c=0;c<nmove;c++){
        unsigned ref;
        tdouble3 mvsimple;
        tmatrix4d mvmatrix;
        if(Motion->GetMov(c,ref,mvsimple,mvmatrix)){//-Movimiento simple
          const unsigned pini=MotionObjBegin[ref]-CaseNfixed,np=MotionObjBegin[ref+1]-MotionObjBegin[ref];
          mvsimple=OrderCode(mvsimple);
          if(Simulate2D)mvsimple.y=0;
          const tfloat3 mvvel=ToTFloat3(mvsimple/TDouble3(stepdt));
          cusph::MoveLinBound(PeriActive,np,pini,mvsimple,mvvel,RidpMoveg,Posxyg,Poszg,Dcellg,Velrhopg,Codeg);
        }
        else{//-Movimiento con matriz
          const unsigned pini=MotionObjBegin[ref]-CaseNfixed,np=MotionObjBegin[ref+1]-MotionObjBegin[ref];
          mvmatrix=OrderCode(mvmatrix);
          //sprintf(Cad,"--->pini:%u np:%u",pini,np); Log->Print(Cad);
          cusph::MoveMatBound(PeriActive,Simulate2D,np,pini,mvmatrix,stepdt,RidpMoveg,Posxyg,Poszg,Dcellg,Velrhopg,Codeg);
        } 
      }
      BoundChanged=true;
    }
  }
  //-Procesa otros modos de motion.
  if(WaveGen){
    if(!nmove)cusph::CalcRidp(PeriActive!=0,Npb,0,CaseNfixed,CaseNfixed+CaseNmoving,Codeg,Idpg,RidpMoveg);
    BoundChanged=true;
    //-Gestion de WaveGen.
    if(WaveGen)for(unsigned c=0;c<WaveGen->GetCount();c++){
      tdouble3 mvsimple;
      tmatrix4d mvmatrix;
      unsigned nparts;
      unsigned idbegin;
      if(WaveGen->GetMotion(c,TimeStep+MotionTimeMod,stepdt,mvsimple,mvmatrix,nparts,idbegin)){//-Movimiento simple
        mvsimple=OrderCode(mvsimple);
        if(Simulate2D)mvsimple.y=0;
        const tfloat3 mvvel=ToTFloat3(mvsimple/TDouble3(stepdt));
        cusph::MoveLinBound(PeriActive,nparts,idbegin-CaseNfixed,mvsimple,mvvel,RidpMoveg,Posxyg,Poszg,Dcellg,Velrhopg,Codeg);
      }
      else{
        mvmatrix=OrderCode(mvmatrix);
        cusph::MoveMatBound(PeriActive,Simulate2D,nparts,idbegin-CaseNfixed,mvmatrix,stepdt,RidpMoveg,Posxyg,Poszg,Dcellg,Velrhopg,Codeg);
      }
    }
  }
  TmgStop(Timers,TMG_SuMotion);
}

//==============================================================================
// Muestra los temporizadores activos.
//==============================================================================
void JSphGpu::ShowTimers(bool onlyfile){
  JLog2::TpMode_Out mode=(onlyfile? JLog2::Out_File: JLog2::Out_ScrFile);
  Log->Print("\n[GPU Timers]",mode);
  if(!SvTimers)Log->Print("none",mode);
  else for(unsigned c=0;c<TimerGetCount();c++)if(TimerIsActive(c))Log->Print(TimerToText(c),mode);
}

//==============================================================================
// Devuelve string con nombres y valores de los timers activos.
//==============================================================================
void JSphGpu::GetTimersInfo(std::string &hinfo,std::string &dinfo)const{
  for(unsigned c=0;c<TimerGetCount();c++)if(TimerIsActive(c)){
    hinfo=hinfo+";"+TimerGetName(c);
    dinfo=dinfo+";"+fun::FloatStr(TimerGetValue(c)/1000.f);
  }
}

//==============================================================================
// Graba fichero vtk con datos de las particulas.
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin,bool idp,bool vel,bool rhop,bool code)
{
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  if(Log->GetMpiRank()>=0)filename=string("p")+fun::IntStr(Log->GetMpiRank())+"_"+filename;
  filename=DirOut+filename;
  //-Reserva memoria basica.
  const unsigned n=pfin-pini;
  ParticlesDataDown(n,pini,code,true,false);
  tfloat3 *pos=new tfloat3[n];
  for(unsigned p=0;p<n;p++)pos[p]=ToTFloat3(AuxPos[p]);
  byte *type=new byte[n];
  for(unsigned p=0;p<n;p++){
    word cod=Code[p];
    byte tp=99;
    if(CODE_GetType(cod)==CODE_TYPE_FIXED)tp=0;
    else if(CODE_GetType(cod)==CODE_TYPE_MOVING)tp=1;
    else if(CODE_GetType(cod)==CODE_TYPE_FLOATING)tp=2;
    else if(CODE_GetType(cod)==CODE_TYPE_FLUID)tp=3;
    if(CODE_GetSpecialValue(cod)==CODE_NORMAL)tp+=0;
    else if(CODE_GetSpecialValue(cod)==CODE_PERIODIC)tp+=10;
    else if(CODE_GetSpecialValue(cod)==CODE_OUTIGNORE)tp+=20;
    else tp+=30;
    type[p]=tp;
  }
  JFormatFiles2::StScalarData fields[8];
  unsigned nfields=0;
  if(idp){  fields[nfields]=JFormatFiles2::DefineField("Id"  ,JFormatFiles2::UInt32  ,1,Idp);     nfields++; }
  if(type){ fields[nfields]=JFormatFiles2::DefineField("Typex",JFormatFiles2::UChar8 ,1,type);    nfields++; }
  if(vel){  fields[nfields]=JFormatFiles2::DefineField("Vel" ,JFormatFiles2::Float32 ,3,AuxVel);  nfields++; }
  if(rhop){ fields[nfields]=JFormatFiles2::DefineField("Rhop",JFormatFiles2::Float32 ,1,AuxRhop); nfields++; }
  if(code){ fields[nfields]=JFormatFiles2::DefineField("Code",JFormatFiles2::UShort16,1,Code);    nfields++; }
  //-Genera fichero.
  JFormatFiles2::SaveVtk(filename,n,pos,nfields,fields);
  //-Libera memoria. 
  delete[] pos;
  delete[] type;
}

//==============================================================================
// Graba fichero vtk con datos de las particulas.
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin,const float3 *posg,const byte *checkg,const unsigned *idpg,const float3 *velg,const float *rhopg){
  //-Reserva memoria basica.
  const unsigned n=pfin-pini;
  tfloat3 *pos=new tfloat3[n];
  byte *check=NULL;
  unsigned *idp=NULL;
  tfloat3 *vel=NULL;
  float *rhop=NULL;
  if(checkg)check=new byte[n];
  if(idpg)idp=new unsigned[n];
  if(velg)vel=new tfloat3[n];
  if(rhopg)rhop=new float[n];
  //-Recupera datos de Gpu.
  cudaMemcpy(pos,posg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(checkg)cudaMemcpy(check,checkg+pini,sizeof(byte)*n,cudaMemcpyDeviceToHost);
  if(idpg)cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  if(velg)cudaMemcpy(vel,velg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(rhopg)cudaMemcpy(rhop,rhopg+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  //-Genera fichero vtk.
  DgSaveVtkParticlesCpu(filename,numfile,0,n,pos,check,idp,vel,rhop);
  //-Libera memoria 
  delete[] pos;
  delete[] check;
  delete[] idp;
  delete[] vel;
  delete[] rhop;
}

//==============================================================================
// Graba fichero csv con datos de las particulas.
//==============================================================================
void JSphGpu::DgSaveCsvParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin,std::string head,const float3 *posg,const unsigned *idpg,const float3 *velg,const float *rhopg,const float *arg,const float3 *aceg,const float3 *vcorrg){
  const char met[]="DgSaveCsvParticlesGpu";
  //-Reserva memoria basica.
  const unsigned n=pfin-pini;
  unsigned *idp=NULL;  if(idpg)idp=new unsigned[n];
  tfloat3 *pos=NULL;   if(posg)pos=new tfloat3[n];
  tfloat3 *vel=NULL;   if(velg)vel=new tfloat3[n];
  float *rhop=NULL;    if(rhopg)rhop=new float[n];
  float *ar=NULL;      if(arg)ar=new float[n];
  tfloat3 *ace=NULL;   if(aceg)ace=new tfloat3[n];
  tfloat3 *vcorr=NULL; if(vcorrg)vcorr=new tfloat3[n];
  //-Recupera datos de Gpu.
  if(idpg)cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  if(posg)cudaMemcpy(pos,posg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(velg)cudaMemcpy(vel,velg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(rhopg)cudaMemcpy(rhop,rhopg+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  if(arg)cudaMemcpy(ar,arg+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  if(aceg)cudaMemcpy(ace,aceg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(vcorrg)cudaMemcpy(vcorr,vcorrg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  CheckCudaError(met,"Failed copying data from GPU.");
  //-Genera fichero csv.
  DgSaveCsvParticlesCpu(filename,numfile,0,n,head,pos,idp,vel,rhop,ar,ace,vcorr);
  //-Libera memoria 
  delete[] idp;
  delete[] pos;
  delete[] vel;
  delete[] rhop;
  delete[] ar;
  delete[] ace;
  delete[] vcorr;
}

//==============================================================================
// Graba fichero csv con datos de las particulas.
//==============================================================================
void JSphGpu::DgSaveCsvParticlesGpu2(std::string filename,int numfile,unsigned pini,unsigned pfin,std::string head,const float3 *posg,const unsigned *idpg,const float3 *velg,const float *rhopg,const float4 *pospresg,const float4 *velrhopg){
  const char met[]="DgSaveCsvParticlesGpu2";
  //-Reserva memoria basica.
  const unsigned n=pfin-pini;
  unsigned *idp=NULL;  if(idpg)idp=new unsigned[n];
  tfloat3 *pos=NULL;   if(posg)pos=new tfloat3[n];
  tfloat3 *vel=NULL;   if(velg)vel=new tfloat3[n];
  float *rhop=NULL;    if(rhopg)rhop=new float[n];
  tfloat4 *pospres=NULL;  if(pospresg)pospres=new tfloat4[n];
  tfloat4 *velrhop=NULL;  if(velrhopg)velrhop=new tfloat4[n];
  //-Recupera datos de Gpu.
  if(idpg)cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  if(posg)cudaMemcpy(pos,posg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(velg)cudaMemcpy(vel,velg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(rhopg)cudaMemcpy(rhop,rhopg+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  if(pospresg)cudaMemcpy(pospres,pospresg+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
  if(velrhopg)cudaMemcpy(velrhop,velrhopg+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
  CheckCudaError(met,"Failed copying data from GPU.");
  //-Genera fichero csv.
  DgSaveCsvParticles2(filename,numfile,0,n,head,pos,idp,vel,rhop,pospres,velrhop);
  //-Libera memoria 
  delete[] idp;
  delete[] pos;
  delete[] vel;
  delete[] rhop;
  delete[] pospres;
  delete[] velrhop;
}


//==============================================================================
// Graba fichero csv con datos de las particulas.
//==============================================================================
void JSphGpu::DgSaveCsvParticles2(std::string filename,int numfile,unsigned pini,unsigned pfin,std::string head,const tfloat3 *pos,const unsigned *idp,const tfloat3 *vel,const float *rhop,const tfloat4 *pospres,const tfloat4 *velrhop){
  const char met[]="DgSaveCsvParticlesCpu";
  int mpirank=Log->GetMpiRank();
  if(mpirank>=0)filename=string("p")+fun::IntStr(mpirank)+"_"+filename;
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  filename=DirOut+filename;
  //-Genera fichero CSV.
  ofstream pf;
  pf.open(filename.c_str());
  if(pf){
    if(!head.empty())pf << head << endl;
    pf << "Num";
    if(idp)pf << ";Idp";
    if(pos)pf << ";PosX;PosY;PosZ";
    if(vel)pf << ";VelX;VelY;VelZ";
    if(rhop)pf << ";Rhop";
    if(pospres)pf << ";Px;Py;Pz;Pp";
    if(velrhop)pf << ";Vx;Vy;Vz;Vr";
    pf << endl;
    const char fmt1[]="%f"; //="%24.16f";
    const char fmt3[]="%f;%f;%f"; //="%24.16f;%24.16f;%24.16f";
    for(unsigned p=pini;p<pfin;p++){
      pf << fun::UintStr(p-pini);
      if(idp)pf << ";" << fun::UintStr(idp[p]);
      if(pos)pf << ";" << fun::Float3Str(pos[p],fmt3);
      if(vel)pf << ";" << fun::Float3Str(vel[p],fmt3);
      if(rhop)pf << ";" << fun::FloatStr(rhop[p],fmt1);
      if(pospres)pf << ";" <<  fun::FloatStr(pospres[p].x,fmt1) << ";" << fun::FloatStr(pospres[p].y,fmt1) << ";" << fun::FloatStr(pospres[p].z,fmt1) << ";" << fun::FloatStr(pospres[p].w,fmt1);
      if(velrhop)pf << ";" <<  fun::FloatStr(velrhop[p].x,fmt1) << ";" << fun::FloatStr(velrhop[p].y,fmt1) << ";" << fun::FloatStr(velrhop[p].z,fmt1) << ";" << fun::FloatStr(velrhop[p].w,fmt1);
      pf << endl;
    }
    if(pf.fail())RunException(met,"Failed writing to file.",filename);
    pf.close();
  }
  else RunException(met,"File could not be opened.",filename);
}


