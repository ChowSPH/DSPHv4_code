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

#include "JSphCpu.h"
#include "JCellDivCpu.h"
#include "JPartFloatBi4.h"
#include "Functions.h"
#include "JSphMotion.h"
#include "JArraysCpu.h"
#include "JSphDtFixed.h"
#include "JWaveGen.h"
#include "JXml.h"
#include "JSaveDt.h"

#include "JFormatFiles2.h"
#include <climits>

#ifdef _WITHOMP
  #include <omp.h>  //Activar tb en Props config -> C/C++ -> Lenguaje -> OpenMp
#else
  #define omp_get_thread_num() 0
  #define omp_get_max_threads() 1
#endif

using namespace std;

//==============================================================================
// Constructor.
//==============================================================================
JSphCpu::JSphCpu(bool withmpi):JSph(true,withmpi){
  ClassName="JSphCpu";
  CellDiv=NULL;
  ArraysCpu=new JArraysCpu;
  InitVars();
  TmcCreation(Timers,false);
}

//==============================================================================
// Destructor.
//==============================================================================
JSphCpu::~JSphCpu(){
  FreeCpuMemoryParticles();
  FreeCpuMemoryFixed();
  delete ArraysCpu;
  TmcDestruction(Timers);
}

//==============================================================================
// Initialization of variables.
//==============================================================================
void JSphCpu::InitVars(){
  RunMode="";
  OmpThreads=1;

  Np=Npb=NpbOk=0;
  NpbPer=NpfPer=0;
  WithFloating=false;

  Idpc=NULL; Codec=NULL; Dcellc=NULL; Posc=NULL; Velrhopc=NULL;
  VelrhopM1c=NULL;                //-Verlet
  PosPrec=NULL; VelrhopPrec=NULL; //-Symplectic
  PsPosc=NULL;                    //-Interaccion Pos-Simple.
  SpsTauc=NULL; SpsGradvelc=NULL; //-Laminar+SPS. 
  Arc=NULL; Acec=NULL; Deltac=NULL;
  ShiftPosc=NULL; ShiftDetectc=NULL; //-Shifting.
  Pressc=NULL;
  RidpMove=NULL; 
  FtRidp=NULL;
  FtoForces=NULL;
  FreeCpuMemoryParticles();
  FreeCpuMemoryFixed();
}

//==============================================================================
// Libera memoria fija en Gpu para moving y floating.
//==============================================================================
void JSphCpu::FreeCpuMemoryFixed(){
  MemCpuFixed=0;
  delete[] RidpMove;  RidpMove=NULL;
  delete[] FtRidp;    FtRidp=NULL;
  delete[] FtoForces; FtoForces=NULL;
}

//==============================================================================
// Allocates memory for arrays with fixed size (motion and floating bodies).
//==============================================================================
void JSphCpu::AllocCpuMemoryFixed(){
  MemCpuFixed=0;
  try{
    //-Allocates memory for moving objects.
    if(CaseNmoving){
      RidpMove=new unsigned[CaseNmoving];  MemCpuFixed+=(sizeof(unsigned)*CaseNmoving);
    }
    //-Allocates memory for floating bodies.
    if(CaseNfloat){
      FtRidp=new unsigned[CaseNfloat];     MemCpuFixed+=(sizeof(unsigned)*CaseNfloat);
      FtoForces=new StFtoForces[FtCount];  MemCpuFixed+=(sizeof(StFtoForces)*FtCount);
    }
  }
  catch(const std::bad_alloc){
    RunException("AllocMemoryFixed","Could not allocate the requested memory.");
  }
}

//==============================================================================
// Libera memoria en Gpu para particulas.
//==============================================================================
void JSphCpu::FreeCpuMemoryParticles(){
  CpuParticlesSize=0;
  MemCpuParticles=0;
  ArraysCpu->Reset();
}

//==============================================================================
// Reserva memoria en Gpu para las particulas. 
//==============================================================================
void JSphCpu::AllocCpuMemoryParticles(unsigned np,float over){
  const char* met="AllocCpuMemoryParticles";
  FreeCpuMemoryParticles();
  //-Calcula numero de particulas para las que se reserva memoria.
  const unsigned np2=(over>0? unsigned(over*np): np);
  CpuParticlesSize=np2;
  //-Calcula cuantos arrays.
  ArraysCpu->SetArraySize(np2);
  ArraysCpu->AddArrayCount(JArraysCpu::SIZE_2B,2);  //-code
  ArraysCpu->AddArrayCount(JArraysCpu::SIZE_4B,5);  //-idp,ar,viscdt,dcell,prrhop
  if(TDeltaSph==DELTA_DynamicExt)ArraysCpu->AddArrayCount(JArraysCpu::SIZE_4B,1);  //-delta
  ArraysCpu->AddArrayCount(JArraysCpu::SIZE_12B,1); //-ace
  ArraysCpu->AddArrayCount(JArraysCpu::SIZE_16B,1); //-velrhop
  ArraysCpu->AddArrayCount(JArraysCpu::SIZE_24B,2); //-pos
  if(Psimple)ArraysCpu->AddArrayCount(JArraysCpu::SIZE_12B,1); //-pspos
  if(TStep==STEP_Verlet){
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_16B,1); //-velrhopm1
  }
  else if(TStep==STEP_Symplectic){
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_24B,1); //-pospre
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_16B,1); //-velrhoppre
  }
  if(TVisco==VISCO_LaminarSPS){     
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_24B,1); //-SpsTau,SpsGradvel
  }
  if(TShifting!=SHIFT_None){
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_12B,1); //-shiftpos
  }
  //-Muestra la memoria reservada.
  MemCpuParticles=ArraysCpu->GetAllocMemoryCpu();
  PrintSizeNp(np2,MemCpuParticles);
}

//==============================================================================
// Resizes space in CPU memory for particles.
//==============================================================================
void JSphCpu::ResizeCpuMemoryParticles(unsigned npnew){
  //-Saves current data from CPU.
  unsigned    *idp       =SaveArrayCpu(Np,Idpc);
  word        *code      =SaveArrayCpu(Np,Codec);
  unsigned    *dcell     =SaveArrayCpu(Np,Dcellc);
  tdouble3    *pos       =SaveArrayCpu(Np,Posc);
  tfloat4     *velrhop   =SaveArrayCpu(Np,Velrhopc);
  tfloat4     *velrhopm1 =SaveArrayCpu(Np,VelrhopM1c);
  tdouble3    *pospre    =SaveArrayCpu(Np,PosPrec);
  tfloat4     *velrhoppre=SaveArrayCpu(Np,VelrhopPrec);
  tsymatrix3f *spstau    =SaveArrayCpu(Np,SpsTauc);
  //-Frees pointers.
  ArraysCpu->Free(Idpc);
  ArraysCpu->Free(Codec);
  ArraysCpu->Free(Dcellc);
  ArraysCpu->Free(Posc);
  ArraysCpu->Free(Velrhopc);
  ArraysCpu->Free(VelrhopM1c);
  ArraysCpu->Free(PosPrec);
  ArraysCpu->Free(VelrhopPrec);
  ArraysCpu->Free(SpsTauc);
  //-Resizes CPU memory allocation.
  const double mbparticle=(double(MemCpuParticles)/(1024*1024))/CpuParticlesSize; //-MB por particula.
  Log->Printf("**JSphCpu: Requesting gpu memory for %u particles: %.1f MB.",npnew,mbparticle*npnew);
  ArraysCpu->SetArraySize(npnew);
  //-Reserve pointers.
  Idpc    =ArraysCpu->ReserveUint();
  Codec   =ArraysCpu->ReserveWord();
  Dcellc  =ArraysCpu->ReserveUint();
  Posc    =ArraysCpu->ReserveDouble3();
  Velrhopc=ArraysCpu->ReserveFloat4();
  if(velrhopm1) VelrhopM1c =ArraysCpu->ReserveFloat4();
  if(pospre)    PosPrec    =ArraysCpu->ReserveDouble3();
  if(velrhoppre)VelrhopPrec=ArraysCpu->ReserveFloat4();
  if(spstau)    SpsTauc    =ArraysCpu->ReserveSymatrix3f();
  //-Restore data in CPU memory.
  RestoreArrayCpu(Np,idp,Idpc);
  RestoreArrayCpu(Np,code,Codec);
  RestoreArrayCpu(Np,dcell,Dcellc);
  RestoreArrayCpu(Np,pos,Posc);
  RestoreArrayCpu(Np,velrhop,Velrhopc);
  RestoreArrayCpu(Np,velrhopm1,VelrhopM1c);
  RestoreArrayCpu(Np,pospre,PosPrec);
  RestoreArrayCpu(Np,velrhoppre,VelrhopPrec);
  RestoreArrayCpu(Np,spstau,SpsTauc);
  //-Updates values.
  CpuParticlesSize=npnew;
  MemCpuParticles=ArraysCpu->GetAllocMemoryCpu();
}

//==============================================================================
// Saves a CPU array in CPU memory. 
//==============================================================================
template<class T> T* JSphCpu::TSaveArrayCpu(unsigned np,const T *datasrc)const{
  T *data=NULL;
  if(datasrc){
    try{
      data=new T[np];
    }
    catch(const std::bad_alloc){
      RunException("TSaveArrayCpu","Could not allocate the requested memory.");
    }
    memcpy(data,datasrc,sizeof(T)*np);
  }
  return(data);
}

//==============================================================================
// Restores a GPU array from CPU memory. 
//==============================================================================
template<class T> void JSphCpu::TRestoreArrayCpu(unsigned np,T *data,T *datanew)const{
  if(data&&datanew)memcpy(datanew,data,sizeof(T)*np);
  delete[] data;
}
//==============================================================================
void JSphCpu::RestoreArrayCpu_Uint(unsigned np,unsigned *data,unsigned *datanew)const{
  if(data&&datanew)memcpy(datanew,data,sizeof(unsigned)*np);
  delete[] data;
}

//==============================================================================
// Arrays para datos basicos de las particulas. 
//==============================================================================
void JSphCpu::ReserveBasicArraysCpu(){
  Idpc=ArraysCpu->ReserveUint();
  Codec=ArraysCpu->ReserveWord();
  Dcellc=ArraysCpu->ReserveUint();
  Posc=ArraysCpu->ReserveDouble3();
  Velrhopc=ArraysCpu->ReserveFloat4();
  if(TStep==STEP_Verlet)VelrhopM1c=ArraysCpu->ReserveFloat4();
  if(TVisco==VISCO_LaminarSPS)SpsTauc=ArraysCpu->ReserveSymatrix3f();
}

//==============================================================================
// Devuelve la memoria reservada en cpu.
//==============================================================================
llong JSphCpu::GetAllocMemoryCpu()const{  
  llong s=JSph::GetAllocMemoryCpu();
  //Reservada en AllocCpuMemoryParticles()
  s+=MemCpuParticles;
  //Reservada en AllocCpuMemoryFixed()
  s+=MemCpuFixed;
  //Reservada en otros objetos
  return(s);
}

//==============================================================================
// Visualiza la memoria reservada
//==============================================================================
void JSphCpu::PrintAllocMemory(llong mcpu)const{
  Log->Printf("Allocated memory in CPU: %lld (%.2f MB)",mcpu,double(mcpu)/(1024*1024));
}

//==============================================================================
// Recupera datos de un rango de particulas y devuelve el numero de particulas que
// sera menor que n si se eliminaron las periodicas.
// - cellorderdecode: Reordena componentes de pos y vel segun CellOrder.
// - onlynormal: Solo se queda con las normales, elimina las particulas periodicas.
//==============================================================================
unsigned JSphCpu::GetParticlesData(unsigned n,unsigned pini,bool cellorderdecode,bool onlynormal
  ,unsigned *idp,tdouble3 *pos,tfloat3 *vel,float *rhop,word *code)
{
  const char met[]="GetParticlesData";
  unsigned num=n;
  //-Copia datos seleccionados.
  if(code)memcpy(code,Codec+pini,sizeof(word)*n);
  if(idp)memcpy(idp,Idpc+pini,sizeof(unsigned)*n);
  if(pos)memcpy(pos,Posc+pini,sizeof(tdouble3)*n);
  if(vel && rhop){
    for(unsigned p=0;p<n;p++){
      tfloat4 vr=Velrhopc[p+pini];
      vel[p]=TFloat3(vr.x,vr.y,vr.z);
      rhop[p]=vr.w;
    }
  }
  else{
    if(vel) for(unsigned p=0;p<n;p++){ tfloat4 vr=Velrhopc[p+pini]; vel[p]=TFloat3(vr.x,vr.y,vr.z); }
    if(rhop)for(unsigned p=0;p<n;p++)rhop[p]=Velrhopc[p+pini].w;
  }
  //-Elimina particulas no normales (periodicas y otras).
  if(onlynormal){
    if(!idp || !pos || !vel || !rhop)RunException(met,"Pointers without data.");
    word *code2=code;
    if(!code2){
      code2=ArraysCpu->ReserveWord();
      memcpy(code2,Codec+pini,sizeof(word)*n);
    }
    unsigned ndel=0;
    for(unsigned p=0;p<n;p++){
      bool normal=(CODE_GetSpecialValue(code2[p])==CODE_NORMAL);
      if(ndel && normal){
        const unsigned pdel=p-ndel;
        idp[pdel]  =idp[p];
        pos[pdel]  =pos[p];
        vel[pdel]  =vel[p];
        rhop[pdel] =rhop[p];
        code2[pdel]=code2[p];
      }
      if(!normal)ndel++;
    }
    num-=ndel;
    if(!code)ArraysCpu->Free(code2);
  }
  //-Reordena componentes en su orden original.
  if(cellorderdecode)DecodeCellOrder(n,pos,vel);
  return(num);
}

//==============================================================================
// Carga la configuracion de ejecucion con OpenMP.
//==============================================================================
void JSphCpu::ConfigOmp(const JCfgRun *cfg){
#ifdef _WITHOMP
  //-Determina numero de threads por host con OpenMP
  if(Cpu && cfg->OmpThreads!=1){
    OmpThreads=cfg->OmpThreads;
    if(OmpThreads<=0)OmpThreads=max(omp_get_num_procs(),1);
    if(OmpThreads>MAXTHREADS_OMP)OmpThreads=MAXTHREADS_OMP;
    omp_set_num_threads(OmpThreads);
    Log->Printf("Threads by host for parallel execution: %d",omp_get_max_threads());
  }
  else{
    OmpThreads=1;
    omp_set_num_threads(OmpThreads);
  }
#else
  OmpThreads=1;
#endif
}

//==============================================================================
// Configura modo de ejecucion en CPU.
//==============================================================================
void JSphCpu::ConfigRunMode(const JCfgRun *cfg,std::string preinfo){
  #ifndef WIN32
    const int len=128; char hname[len];
    gethostname(hname,len);
    if(!preinfo.empty())preinfo=preinfo+", ";
    preinfo=preinfo+"HostName:"+hname;
  #endif
  Hardware="Cpu";
  if(OmpThreads==1)RunMode="Single core";
  else RunMode=string("OpenMP(Threads:")+fun::IntStr(OmpThreads)+")";
  if(!preinfo.empty())RunMode=preinfo+", "+RunMode;
  if(Stable)RunMode=string("Stable, ")+RunMode;
  if(Psimple)RunMode=string("Pos-Simple, ")+RunMode;
  else RunMode=string("Pos-Double, ")+RunMode;
  Log->Print(fun::VarStr("RunMode",RunMode));
}

//==============================================================================
// Inicializa vectores y variables para la ejecucion.
//==============================================================================
void JSphCpu::InitRun(){
  const char met[]="InitRun";
  WithFloating=(CaseNfloat>0);
  if(TStep==STEP_Verlet){
    memcpy(VelrhopM1c,Velrhopc,sizeof(tfloat4)*Np);
    VerletStep=0;
  }
  else if(TStep==STEP_Symplectic)DtPre=DtIni;
  if(TVisco==VISCO_LaminarSPS)memset(SpsTauc,0,sizeof(tsymatrix3f)*Np);
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
}

//==============================================================================
// Prepara variables para interaccion "INTER_Forces" o "INTER_ForcesCorr".
//==============================================================================
void JSphCpu::PreInteractionVars_Forces(TpInter tinter,unsigned ini,unsigned np,unsigned npb){
  //-Inicializa arrays.
  const unsigned npf=np-npb;
  if(ini)RunException("PreInteractionVars_Forces","initialization of arrays is invalid.");
  memset(Arc+ini,0,sizeof(float)*np);                                    //Arc[]=0
  if(Deltac)memset(Deltac+ini,0,sizeof(float)*np);                       //Deltac[]=0
  if(ShiftPosc)memset(ShiftPosc+ini,0,sizeof(tfloat3)*np);               //ShiftPosc[]=0
  if(ShiftDetectc)memset(ShiftDetectc+ini,0,sizeof(float)*np);           //ShiftDetectc[]=0
  memset(Acec+ini,0,sizeof(tfloat3)*npb);                                //Acec[]=(0,0,0) para bound
  for(unsigned p=npb;p<np;p++)Acec[p]=Gravity;                           //Acec[]=Gravity para fluid
  if(SpsGradvelc)memset(SpsGradvelc+ini+npb,0,sizeof(tsymatrix3f)*npf);  //SpsGradvelc[]=(0,0,0,0,0,0).

  //-Prepara datos derivados de rhop para interaccion.
  const int n=int(np);
  #ifdef _WITHOMP
    #pragma omp parallel for schedule (static) if(n>LIMIT_PREINTERACTION_OMP)
  #endif
  for(int p=0;p<n;p++){
    const float rhop=Velrhopc[p].w,rhop_r0=rhop*OVERRHOPCERO;
    Pressc[p]=CteB*(pow(rhop_r0,Gamma)-1.0f);
  }
}

//==============================================================================
// Prepara variables para interaccion "INTER_Forces" o "INTER_ForcesCorr".
//==============================================================================
void JSphCpu::PreInteraction_Forces(TpInter tinter){
  TmcStart(Timers,TMC_CfPreForces);
  //-Asigna memoria.
  Arc=ArraysCpu->ReserveFloat();
  Acec=ArraysCpu->ReserveFloat3();
  if(TDeltaSph==DELTA_DynamicExt)Deltac=ArraysCpu->ReserveFloat();
  if(TShifting!=SHIFT_None){
    ShiftPosc=ArraysCpu->ReserveFloat3();
    if(ShiftTFS)ShiftDetectc=ArraysCpu->ReserveFloat();
  }   
  Pressc=ArraysCpu->ReserveFloat();
  if(TVisco==VISCO_LaminarSPS)SpsGradvelc=ArraysCpu->ReserveSymatrix3f();

  //-Prepara datos para interaccion Pos-Simple.
  if(Psimple){
    PsPosc=ArraysCpu->ReserveFloat3();
    const int np=int(Np);
    #ifdef _WITHOMP
      #pragma omp parallel for schedule (static) if(np>LIMIT_PREINTERACTION_OMP)
    #endif
    for(int p=0;p<np;p++){ PsPosc[p]=ToTFloat3(Posc[p]); }
  }
  //-Inicializa arrays.
  PreInteractionVars_Forces(tinter,0,Np,Npb);

    //-Calcula VelMax: Se incluyen las particulas floatings y no afecta el uso de condiciones periodicas.
#ifdef DT_ALLPARTICLES
    const unsigned pini=0;
  #else
    const unsigned pini=Npb;
  #endif
  float velmax=0;
  for(unsigned p=pini;p<Np;p++){
    const tfloat4 v=Velrhopc[p];
    const float v2=v.x*v.x+v.y*v.y+v.z*v.z;
    velmax=max(velmax,v2);
  }
  VelMax=sqrt(velmax);
  ViscDtMax=0;
  TmcStop(Timers,TMC_CfPreForces);
}

//==============================================================================
// Libera memoria asignada de ArraysGpu.
//==============================================================================
void JSphCpu::PosInteraction_Forces(){
  //-Libera memoria asignada en PreInteraction_Forces().
  ArraysCpu->Free(Arc);          Arc=NULL;
  ArraysCpu->Free(Acec);         Acec=NULL;
  ArraysCpu->Free(Deltac);       Deltac=NULL;
  ArraysCpu->Free(ShiftPosc);    ShiftPosc=NULL;
  ArraysCpu->Free(ShiftDetectc); ShiftDetectc=NULL;
  ArraysCpu->Free(Pressc);       Pressc=NULL;
  ArraysCpu->Free(PsPosc);       PsPosc=NULL;
  ArraysCpu->Free(SpsGradvelc);  SpsGradvelc=NULL;
}

//==============================================================================
// Realiza interaccion entre particulas para Ren correction. Bound-Fluid/Float
//==============================================================================
template<bool psimple,TpFtMode ftmode> void JSphCpu::InteractionRenBound
  (unsigned n,unsigned pinit,tint4 nc,int hdiv,unsigned cellinitial
  ,const unsigned *beginendcell,tint3 cellzero,const unsigned *dcell
  ,const tdouble3 *pos,const tfloat3 *pspos,const tfloat4 *velrhop,const word *code,const unsigned *idp
  ,const float *press,float *presskf)
{
  //-Inicia ejecucion con OpenMP.
  const int pfin=int(pinit+n);
  #ifdef _WITHOMP
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=int(pinit);p1<pfin;p1++){
    float pkfap1=0,pkfbp1=0;
    const tfloat3 psposp1=(psimple? pspos[p1]: TFloat3(0));
    const tdouble3 posp1=(psimple? TDouble3(0): pos[p1]);
    const float rhopp1=velrhop[p1].w;
    const float pressp1=press[p1];

    //-Obtiene limites de interaccion
    unsigned rcell=dcell[p1];
    const int cx=PC__Cellx(DomCellCode,rcell)-cellzero.x;
    const int cy=PC__Celly(DomCellCode,rcell)-cellzero.y;
    const int cz=PC__Cellz(DomCellCode,rcell)-cellzero.z;
    //-Codigo para hdiv 1 o 2 pero no cero.
    const int cxini=cx-min(cx,hdiv);
    const int cxfin=cx+min(nc.x-cx-1,hdiv)+1;
    const int yini=cy-min(cy,hdiv);
    const int yfin=cy+min(nc.y-cy-1,hdiv)+1;
    const int zini=cz-min(cz,hdiv);
    const int zfin=cz+min(nc.z-cz-1,hdiv)+1;

    //-Busqueda de vecinos en celdas adyacentes.
    for(int z=zini;z<zfin;z++){
      const int zmod=(nc.w)*z+cellinitial; //-Le suma donde empiezan las celdas de fluido.
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+nc.x*y;
        const unsigned pini=beginendcell[cxini+ymod];
        const unsigned pfin=beginendcell[cxfin+ymod];

        //-Interaccion de Bound con varias Fluid/Float.
        //----------------------------------------------
        for(unsigned p2=pini;p2<pfin;p2++){
          const float drx=(psimple? psposp1.x-pspos[p2].x: float(posp1.x-pos[p2].x));
          const float dry=(psimple? psposp1.y-pspos[p2].y: float(posp1.y-pos[p2].y));
          const float drz=(psimple? psposp1.z-pspos[p2].z: float(posp1.z-pos[p2].z));
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=Fourh2 && rr2>=1e-18f){
            float wab;
            {//-Wendland kernel
              const float qq=sqrt(rr2)/H;
              const float wqq=2.f*qq+1.f;
              const float wqq1=1.f-0.5f*qq;
              const float wqq2=wqq1*wqq1;
              wab=Awen*wqq*wqq2*wqq2;
            }

            //===== Obtiene datos de particula p2 ===== 
            const float rhopp2=velrhop[p2].w;
            float massp2=MassFluid; //-Contiene masa de particula por defecto fluid.
            if(USE_FLOATING){
              bool ftp2=(CODE_GetType(code[p2])==CODE_TYPE_FLOATING);
              if(ftp2)massp2=FtObjs[CODE_GetTypeValue(code[p2])].massp;
            }

            //-Acumula numerador y denominador para calculo de Pkf.
            const float pkf=(massp2/rhopp2)*wab;
            pkfap1+=pkf*(press[p2]+rhopp2*Gravity.z*drz); //<--La gravedad debe aplicarse de forma general...
            pkfbp1+=pkf;
          }
        }
      }
    }
    //-Almacena resultados.
    presskf[p1]=(pkfbp1!=0? pkfap1/pkfbp1: 0); //<--Se deberia controlar cuando se aplica en funci�n de las particulas de fluido que vea...
  }
}

//==============================================================================
// Calcula presion suavizada de fluido en contorno para Ren Correction. Bound-Fluid/Float
//==============================================================================
void JSphCpu::Interaction_Ren(unsigned np,unsigned npb,unsigned npbok
  ,tuint3 ncells,const unsigned *begincell,tuint3 cellmin,const unsigned *dcell
  ,const tdouble3 *pos,const tfloat3 *pspos,const tfloat4 *velrhop,const unsigned *idp,const word *code
  ,const float *press,float *presskf)
{
  const tint4 nc=TInt4(int(ncells.x),int(ncells.y),int(ncells.z),int(ncells.x*ncells.y));
  const tint3 cellzero=TInt3(cellmin.x,cellmin.y,cellmin.z);
  const unsigned cellfluid=nc.w*nc.z+1;
  const int hdiv=(CellMode==CELLMODE_H? 2: 1);
  if(!WithFloating){ const TpFtMode ftmode=FTMODE_None;
    if(Psimple)InteractionRenBound<true ,ftmode> (npbok,0,nc,hdiv,cellfluid,begincell,cellzero,dcell,pos,pspos,velrhop,code,idp,press,presskf);
    else       InteractionRenBound<false,ftmode> (npbok,0,nc,hdiv,cellfluid,begincell,cellzero,dcell,pos,pspos,velrhop,code,idp,press,presskf);
  }
  else{              const TpFtMode ftmode=FTMODE_Sph;  
    if(Psimple)InteractionRenBound<true ,ftmode> (npbok,0,nc,hdiv,cellfluid,begincell,cellzero,dcell,pos,pspos,velrhop,code,idp,press,presskf);
    else       InteractionRenBound<false,ftmode> (npbok,0,nc,hdiv,cellfluid,begincell,cellzero,dcell,pos,pspos,velrhop,code,idp,press,presskf);
  }
}

//==============================================================================
// Calcula nuevo valor de presion y desidad aplicando Ren correction.
//==============================================================================
void JSphCpu::ComputeRenPress(unsigned npbok,float beta,const float *presskf,tfloat4 *velrhop,float *press)const{
  const int n=int(npbok);
  #ifdef _WITHOMP
    #pragma omp parallel for schedule (static) if(n>LIMIT_PREINTERACTION_OMP)
  #endif
  for(int p=0;p<n;p++){
    const float pressc=press[p]=beta*presskf[p]+(1.f-beta)*press[p];
    Velrhopc[p].w=RhopZero*pow(pressc/CteB+1.f,1.f/Gamma);
  }
}



//==============================================================================
// Realiza interaccion entre particulas. Bound-Fluid/Float
//==============================================================================
template<bool psimple,TpFtMode ftmode> void JSphCpu::InteractionForcesBound
  (unsigned n,unsigned pinit,tint4 nc,int hdiv,unsigned cellinitial
  ,const unsigned *beginendcell,tint3 cellzero,const unsigned *dcell
  ,const tdouble3 *pos,const tfloat3 *pspos,const tfloat4 *velrhop,const word *code,const unsigned *idp
  ,float &viscdt,float *ar)
{
  //-Inicializa viscth para calcular visdt maximo con OpenMP.
  float viscth[MAXTHREADS_OMP*STRIDE_OMP];
  for(int th=0;th<OmpThreads;th++)viscth[th*STRIDE_OMP]=0;
  //-Inicia ejecucion con OpenMP.
  const int pfin=int(pinit+n);
  #ifdef _WITHOMP
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=int(pinit);p1<pfin;p1++){
    float visc=0,arp1=0;
    const tfloat3 velp1=TFloat3(velrhop[p1].x,velrhop[p1].y,velrhop[p1].z);
    const tfloat3 psposp1=(psimple? pspos[p1]: TFloat3(0));
    const tdouble3 posp1=(psimple? TDouble3(0): pos[p1]);

    //-Obtiene limites de interaccion
    unsigned rcell=dcell[p1];
    const int cx=PC__Cellx(DomCellCode,rcell)-cellzero.x;
    const int cy=PC__Celly(DomCellCode,rcell)-cellzero.y;
    const int cz=PC__Cellz(DomCellCode,rcell)-cellzero.z;
    //-Codigo para hdiv 1 o 2 pero no cero.
    const int cxini=cx-min(cx,hdiv);
    const int cxfin=cx+min(nc.x-cx-1,hdiv)+1;
    const int yini=cy-min(cy,hdiv);
    const int yfin=cy+min(nc.y-cy-1,hdiv)+1;
    const int zini=cz-min(cz,hdiv);
    const int zfin=cz+min(nc.z-cz-1,hdiv)+1;

    //-Busqueda de vecinos en celdas adyacentes.
    for(int z=zini;z<zfin;z++){
      const int zmod=(nc.w)*z+cellinitial; //-Le suma donde empiezan las celdas de fluido.
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+nc.x*y;
        const unsigned pini=beginendcell[cxini+ymod];
        const unsigned pfin=beginendcell[cxfin+ymod];

        //-Interaccion de Bound con varias Fluid/Float.
        //----------------------------------------------
        for(unsigned p2=pini;p2<pfin;p2++){
          const float drx=(psimple? psposp1.x-pspos[p2].x: float(posp1.x-pos[p2].x));
          const float dry=(psimple? psposp1.y-pspos[p2].y: float(posp1.y-pos[p2].y));
          const float drz=(psimple? psposp1.z-pspos[p2].z: float(posp1.z-pos[p2].z));
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=Fourh2 && rr2>=1e-18f){

            float frx,fry,frz;
            {//===== Kernel =====
              const float rad=sqrt(rr2);
              const float qq=rad/H;
              //-Wendland kernel
              const float wqq1=1.f-0.5f*qq;
              const float wqq2=wqq1*wqq1;
              const float fac=Bwen*qq*wqq2*wqq1/rad;
              frx=fac*drx; fry=fac*dry; frz=fac*drz;
            }

            //===== Obtiene masa de particula p2 ===== 
            float massp2=MassFluid; //-Contiene masa de particula por defecto fluid.
            bool compute=true;      //-Se desactiva cuando se usa DEM y es bound-float.
            if(USE_FLOATING){
              bool ftp2=(CODE_GetType(code[p2])==CODE_TYPE_FLOATING);
              if(ftp2)massp2=FtObjs[CODE_GetTypeValue(code[p2])].massp;
              compute=!(USE_DEM && ftp2); //-Se desactiva cuando se usa DEM y es bound-float.
            }

            if(compute){
              //-Density derivative
              const float dvx=velp1.x-velrhop[p2].x, dvy=velp1.y-velrhop[p2].y, dvz=velp1.z-velrhop[p2].z;
              if(compute)arp1+=massp2*(dvx*frx+dvy*fry+dvz*frz);

              {//===== Viscosity ===== 
                const float dot=drx*dvx + dry*dvy + drz*dvz;
                const float dot_rr2=dot/(rr2+Eta2);
                visc=max(dot_rr2,visc);
              }
            }
          }
        }
      }
    }
    //-Almacena resultados.
    if(arp1||visc){
      ar[p1]+=arp1;
      const int th=omp_get_thread_num();
      if(visc>viscth[th*STRIDE_OMP])viscth[th*STRIDE_OMP]=visc;
    }
  }
  //-Guarda en viscdt el valor maximo.
  for(int th=0;th<OmpThreads;th++)if(viscdt<viscth[th*STRIDE_OMP])viscdt=viscth[th*STRIDE_OMP];
}

//==============================================================================
// Realiza interaccion entre particulas. Fluid/Float-Fluid/Float or Fluid/Float-Bound
//==============================================================================
template<bool psimple,TpFtMode ftmode,bool lamsps,TpDeltaSph tdelta,bool shift> void JSphCpu::InteractionForcesFluid
  (unsigned n,unsigned pinit,tint4 nc,int hdiv,unsigned cellinitial,float visco
  ,const unsigned *beginendcell,tint3 cellzero,const unsigned *dcell
  ,const tsymatrix3f* tau,tsymatrix3f* gradvel
  ,const tdouble3 *pos,const tfloat3 *pspos,const tfloat4 *velrhop,const word *code,const unsigned *idp
  ,const float *press 
  ,float &viscdt,float *ar,tfloat3 *ace,float *delta
  ,TpShifting tshifting,tfloat3 *shiftpos,float *shiftdetect)const
{
  const bool boundp2=(!cellinitial); //-Interaccion con Bound.
  //-Inicializa viscth para calcular visdt maximo con OpenMP.
  float viscth[MAXTHREADS_OMP*STRIDE_OMP];
  for(int th=0;th<OmpThreads;th++)viscth[th*STRIDE_OMP]=0;
  //-Inicia ejecucion con OpenMP.
  const int pfin=int(pinit+n);
  #ifdef _WITHOMP
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=int(pinit);p1<pfin;p1++){
    float visc=0,arp1=0,deltap1=0;
//    bool withbound=false;
    tfloat3 acep1=TFloat3(0);
    tsymatrix3f gradvelp1={0,0,0,0,0,0};
    tfloat3 shiftposp1=TFloat3(0);
    float shiftdetectp1=0;

    //-Obtiene datos de particula p1 en caso de existir floatings.
    bool ftp1=false;     //-Indica si es floating.
    float ftmassp1=1.f;  //-Contiene masa de particula floating o 1.0f si es fluid.
    if(USE_FLOATING){
      ftp1=(CODE_GetType(code[p1])==CODE_TYPE_FLOATING);
      if(ftp1)ftmassp1=FtObjs[CODE_GetTypeValue(code[p1])].massp;
      if(ftp1 && (tdelta==DELTA_Dynamic || tdelta==DELTA_DynamicExt))deltap1=FLT_MAX;
      if(ftp1 && shift)shiftposp1.x=FLT_MAX;  //-Para floatings no se calcula shifting.
    }

    //-Obtiene datos de particula p1.
    const tfloat3 velp1=TFloat3(velrhop[p1].x,velrhop[p1].y,velrhop[p1].z);
    const float rhopp1=velrhop[p1].w;
    const tfloat3 psposp1=(psimple? pspos[p1]: TFloat3(0));
    const tdouble3 posp1=(psimple? TDouble3(0): pos[p1]);
    const float pressp1=press[p1];
    const tsymatrix3f taup1=(lamsps? tau[p1]: gradvelp1);

    //-Obtiene limites de interaccion
    unsigned rcell=dcell[p1];
    const int cx=PC__Cellx(DomCellCode,rcell)-cellzero.x;
    const int cy=PC__Celly(DomCellCode,rcell)-cellzero.y;
    const int cz=PC__Cellz(DomCellCode,rcell)-cellzero.z;

    //-Codigo para hdiv 1 o 2 pero no cero.
    const int cxini=cx-min(cx,hdiv);
    const int cxfin=cx+min(nc.x-cx-1,hdiv)+1;
    const int yini=cy-min(cy,hdiv);
    const int yfin=cy+min(nc.y-cy-1,hdiv)+1;
    const int zini=cz-min(cz,hdiv);
    const int zfin=cz+min(nc.z-cz-1,hdiv)+1;

    //-Busqueda de vecinos en celdas adyacentes.
    for(int z=zini;z<zfin;z++){
      const int zmod=(nc.w)*z+cellinitial; //-Le suma donde empiezan las celdas de fluido o bound.
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+nc.x*y;
        const unsigned pini=beginendcell[cxini+ymod];
        const unsigned pfin=beginendcell[cxfin+ymod];

        //-Interaccion de Fluid con varias Fluid o Bound.
        //------------------------------------------------
        for(unsigned p2=pini;p2<pfin;p2++){
          const float drx=(psimple? psposp1.x-pspos[p2].x: float(posp1.x-pos[p2].x));
          const float dry=(psimple? psposp1.y-pspos[p2].y: float(posp1.y-pos[p2].y));
          const float drz=(psimple? psposp1.z-pspos[p2].z: float(posp1.z-pos[p2].z));
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=Fourh2 && rr2>=1e-18f){
//            //if(boundp2 && CODE_GetType(code[p2])==CODE_TYPE_FIXED)withbound=true;
//            if(boundp2)withbound=true;

            float frx,fry,frz;
            {//===== Kernel =====
              const float rad=sqrt(rr2);
              const float qq=rad/H;
              //-Wendland kernel
              const float wqq1=1.f-0.5f*qq;
              const float wqq2=wqq1*wqq1;
              const float fac=Bwen*qq*wqq2*wqq1/rad;
              frx=fac*drx; fry=fac*dry; frz=fac*drz;
            }

            //===== Obtiene masa de particula p2 ===== 
            float massp2=(boundp2? MassBound: MassFluid); //-Contiene masa de particula segun sea bound o fluid.
            bool ftp2=false;    //-Indica si es floating.
            bool compute=true;  //-Se desactiva cuando se usa DEM y es float-float o float-bound.
            if(USE_FLOATING){
              ftp2=(CODE_GetType(code[p2])==CODE_TYPE_FLOATING);
              if(ftp2)massp2=FtObjs[CODE_GetTypeValue(code[p2])].massp;
              if(ftp2 && (tdelta==DELTA_Dynamic || tdelta==DELTA_DynamicExt))deltap1=FLT_MAX;
              if(ftp2 && shift && tshifting==SHIFT_NoBound)shiftposp1.x=FLT_MAX; //-Con floatings anula shifting.
              compute=!(USE_DEM && ftp1 && (boundp2 || ftp2)); //-Se desactiva cuando se usa DEM y es float-float o float-bound.
            }

            //===== Aceleration ===== 
            if(compute){
              const float prs=(pressp1+press[p2])/(rhopp1*velrhop[p2].w);
              const float p_vpm=-prs*massp2*ftmassp1;
              acep1.x+=p_vpm*frx; acep1.y+=p_vpm*fry; acep1.z+=p_vpm*frz;
            }

            //-Density derivative
            const float dvx=velp1.x-velrhop[p2].x, dvy=velp1.y-velrhop[p2].y, dvz=velp1.z-velrhop[p2].z;
            if(compute)arp1+=massp2*(dvx*frx+dvy*fry+dvz*frz);

            const float cbar=(float)Cs0;
            //-Density derivative (DeltaSPH Molteni)
            if((tdelta==DELTA_Dynamic || tdelta==DELTA_DynamicExt) && deltap1!=FLT_MAX){
              const float rhop1over2=rhopp1/velrhop[p2].w;
              const float visc_densi=Delta2H*cbar*(rhop1over2-1.f)/(rr2+Eta2);
              const float dot3=(drx*frx+dry*fry+drz*frz);
              const float delta=visc_densi*dot3*massp2;
              deltap1=(boundp2? FLT_MAX: deltap1+delta);
            }

            //-Shifting correction
            if(shift && shiftposp1.x!=FLT_MAX){
              const float massrhop=massp2/velrhop[p2].w;
              const bool noshift=(boundp2 && (tshifting==SHIFT_NoBound || (tshifting==SHIFT_NoFixed && CODE_GetType(code[p2])==CODE_TYPE_FIXED)));
              shiftposp1.x=(noshift? FLT_MAX: shiftposp1.x+massrhop*frx); //-Con boundary anula shifting.
              shiftposp1.y+=massrhop*fry;
              shiftposp1.z+=massrhop*frz;
              shiftdetectp1-=massrhop*(drx*frx+dry*fry+drz*frz);
            }

            //===== Viscosity ===== 
            if(compute){
              const float dot=drx*dvx + dry*dvy + drz*dvz;
              const float dot_rr2=dot/(rr2+Eta2);
              visc=max(dot_rr2,visc);
              if(!lamsps){//-Artificial viscosity 
                if(dot<0){
                  const float amubar=H*dot_rr2;  //amubar=CTE.h*dot/(rr2+CTE.eta2);
                  const float robar=(rhopp1+velrhop[p2].w)*0.5f;
                  const float pi_visc=(-visco*cbar*amubar/robar)*massp2*ftmassp1;
                  acep1.x-=pi_visc*frx; acep1.y-=pi_visc*fry; acep1.z-=pi_visc*frz;
                }
              }
              else{//-Laminar+SPS viscosity 
                {//-Laminar contribution.
                  const float robar2=(rhopp1+velrhop[p2].w);
                  const float temp=4.f*visco/((rr2+Eta2)*robar2);  //-Simplificacion de temp=2.0f*visco/((rr2+CTE.eta2)*robar); robar=(rhopp1+velrhop2.w)*0.5f;
                  const float vtemp=massp2*temp*(drx*frx+dry*fry+drz*frz);  
                  acep1.x+=vtemp*dvx; acep1.y+=vtemp*dvy; acep1.z+=vtemp*dvz;
                }
                //-SPS turbulence model.
                float tau_xx=taup1.xx,tau_xy=taup1.xy,tau_xz=taup1.xz; //-taup1 siempre es cero cuando p1 no es fluid.
                float tau_yy=taup1.yy,tau_yz=taup1.yz,tau_zz=taup1.zz;
                if(!boundp2 && !ftp2){//-Cuando p2 es fluido.  
                  tau_xx+=tau[p2].xx; tau_xy+=tau[p2].xy; tau_xz+=tau[p2].xz;
                  tau_yy+=tau[p2].yy; tau_yz+=tau[p2].yz; tau_zz+=tau[p2].zz;
                }
                acep1.x+=massp2*ftmassp1*(tau_xx*frx+tau_xy*fry+tau_xz*frz);
                acep1.y+=massp2*ftmassp1*(tau_xy*frx+tau_yy*fry+tau_yz*frz);
                acep1.z+=massp2*ftmassp1*(tau_xz*frx+tau_yz*fry+tau_zz*frz);
                //-Velocity gradients.
                if(!ftp1){//-Cuando p1 es fluido. 
                  const float volp2=-massp2/velrhop[p2].w;
                  float dv=dvx*volp2; gradvelp1.xx+=dv*frx; gradvelp1.xy+=dv*fry; gradvelp1.xz+=dv*frz;
                        dv=dvy*volp2; gradvelp1.xy+=dv*frx; gradvelp1.yy+=dv*fry; gradvelp1.yz+=dv*frz;
                        dv=dvz*volp2; gradvelp1.xz+=dv*frx; gradvelp1.yz+=dv*fry; gradvelp1.zz+=dv*frz;
                  // to compute tau terms we assume that gradvel.xy=gradvel.dudy+gradvel.dvdx, gradvel.xz=gradvel.dudz+gradvel.dwdx, gradvel.yz=gradvel.dvdz+gradvel.dwdy
                  // so only 6 elements are needed instead of 3x3.
                }
              }
            }
          }
        }
      }
    }
    //-Almacena resultados.
    if(shift||arp1||acep1.x||acep1.y||acep1.z||visc){
      //if(Idpc[p1]==3000)Log->Printf("In-----> p1:%u    px:%f",p1,shiftposp1.x);
      if(tdelta==DELTA_Dynamic&&deltap1!=FLT_MAX)arp1+=deltap1;
      if(tdelta==DELTA_DynamicExt)delta[p1]=(delta[p1]==FLT_MAX || deltap1==FLT_MAX? FLT_MAX: delta[p1]+deltap1);
      ar[p1]+=arp1;
      ace[p1]=ace[p1]+acep1;
      const int th=omp_get_thread_num();
      if(visc>viscth[th*STRIDE_OMP])viscth[th*STRIDE_OMP]=visc;
      if(lamsps){
        gradvel[p1].xx+=gradvelp1.xx;
        gradvel[p1].xy+=gradvelp1.xy;
        gradvel[p1].xz+=gradvelp1.xz;
        gradvel[p1].yy+=gradvelp1.yy;
        gradvel[p1].yz+=gradvelp1.yz;
        gradvel[p1].zz+=gradvelp1.zz;
      }
      if(shift && shiftpos[p1].x!=FLT_MAX){
        shiftpos[p1]=(shiftposp1.x==FLT_MAX? TFloat3(FLT_MAX,0,0): shiftpos[p1]+shiftposp1);
        if(shiftdetect)shiftdetect[p1]+=shiftdetectp1;
      }
    }
  }
  //-Guarda en viscdt el valor maximo.
  for(int th=0;th<OmpThreads;th++)if(viscdt<viscth[th*STRIDE_OMP])viscdt=viscth[th*STRIDE_OMP];
}

//==============================================================================
// Realiza interaccion DEM entre particulas Floating-Bound & Floating-Floating //(DEM)
//==============================================================================
template<bool psimple> void JSphCpu::InteractionForcesDEM
  (unsigned nfloat,tint4 nc,int hdiv,unsigned cellfluid
  ,const unsigned *beginendcell,tint3 cellzero,const unsigned *dcell
  ,const unsigned *ftridp,const StDemData* demobjs
  ,const tdouble3 *pos,const tfloat3 *pspos,const tfloat4 *velrhop,const word *code,const unsigned *idp
  ,float &viscdt,tfloat3 *ace)const
{
  //-Inicializa demdtth para calcular demdt maximo con OpenMP.
  float demdtth[MAXTHREADS_OMP*STRIDE_OMP];
  for(int th=0;th<OmpThreads;th++)demdtth[th*STRIDE_OMP]=-FLT_MAX;
  //-Inicia ejecucion con OpenMP.
  const int nft=int(nfloat);
  #ifdef _WITHOMP
    #pragma omp parallel for schedule (guided)
  #endif
  for(int cf=0;cf<nft;cf++){
    const unsigned p1=ftridp[cf];
    if(p1!=UINT_MAX){
      float demdtp1=0;
      tfloat3 acep1=TFloat3(0);

      //-Obtiene datos de particula p1.
      const tfloat3 psposp1=(psimple? pspos[p1]: TFloat3(0));
      const tdouble3 posp1=(psimple? TDouble3(0): pos[p1]);
      const word tavp1=CODE_GetTypeAndValue(code[p1]);
      const float masstotp1=demobjs[tavp1].mass;
      const float taup1=demobjs[tavp1].tau;
      const float kfricp1=demobjs[tavp1].kfric;
      const float restitup1=demobjs[tavp1].restitu;

      //-Obtiene limites de interaccion
      unsigned rcell=dcell[p1];
      const int cx=PC__Cellx(DomCellCode,rcell)-cellzero.x;
      const int cy=PC__Celly(DomCellCode,rcell)-cellzero.y;
      const int cz=PC__Cellz(DomCellCode,rcell)-cellzero.z;

      //-Codigo para hdiv 1 o 2 pero no cero.
      const int cxini=cx-min(cx,hdiv);
      const int cxfin=cx+min(nc.x-cx-1,hdiv)+1;
      const int yini=cy-min(cy,hdiv);
      const int yfin=cy+min(nc.y-cy-1,hdiv)+1;
      const int zini=cz-min(cz,hdiv);
      const int zfin=cz+min(nc.z-cz-1,hdiv)+1;

      //-Busqueda de vecinos en celdas adyacentes (primero bound y despues fluid+floating).
      for(unsigned cellinitial=0;cellinitial<=cellfluid;cellinitial+=cellfluid){
        for(int z=zini;z<zfin;z++){
          const int zmod=(nc.w)*z+cellinitial; //-Le suma donde empiezan las celdas de fluido o bound.
          for(int y=yini;y<yfin;y++){
            int ymod=zmod+nc.x*y;
            const unsigned pini=beginendcell[cxini+ymod];
            const unsigned pfin=beginendcell[cxfin+ymod];

            //-Interaccion de Floating con varias Fluid o Bound.
            //------------------------------------------------
            for(unsigned p2=pini;p2<pfin;p2++)if(CODE_GetType(code[p2])!=CODE_TYPE_FLUID && tavp1!=CODE_GetTypeAndValue(code[p2])){
              const float drx=(psimple? psposp1.x-pspos[p2].x: float(posp1.x-pos[p2].x));
              const float dry=(psimple? psposp1.y-pspos[p2].y: float(posp1.y-pos[p2].y));
              const float drz=(psimple? psposp1.z-pspos[p2].z: float(posp1.z-pos[p2].z));
              const float rr2=drx*drx+dry*dry+drz*drz;
              const float rad=sqrt(rr2);

              //-Calcula valor maximo de demdt.
              const word tavp2=CODE_GetTypeAndValue(code[p2]);
              const float masstotp2=demobjs[tavp2].mass;
              const float taup2=demobjs[tavp2].tau;
              const float kfricp2=demobjs[tavp2].kfric;
              const float restitup2=demobjs[tavp2].restitu;
                    //const StDemData *demp2=demobjs+CODE_GetTypeAndValue(code[p2]);

                    const float nu_mass=(!cellinitial? masstotp1/2: masstotp1*masstotp2/(masstotp1+masstotp2)); //-Con boundary toma la propia masa del floating 1.
              const float kn=4/(3*(taup1+taup2))*sqrt(float(Dp)/4); //generalized rigidity - Lemieux 2008
                  const float demvisc=float(PI)/(sqrt( kn/nu_mass ))*40.f;              
              if(demdtp1<demvisc)demdtp1=demvisc;

              const float over_lap=1.0f*float(Dp)-rad; //-(ri+rj)-|dij|
                  if(over_lap>0.0f){ //-Contact
                  const float dvx=velrhop[p1].x-velrhop[p2].x, dvy=velrhop[p1].y-velrhop[p2].y, dvz=velrhop[p1].z-velrhop[p2].z; //vji
                  const float nx=drx/rad, ny=dry/rad, nz=drz/rad; //normal_ji               
                  const float vn=dvx*nx+dvy*ny+dvz*nz; //vji.nji      
                      //normal
                const float eij=(restitup1+restitup2)/2;
                      const float gn=-(2.0f*log(eij)*sqrt(nu_mass*kn))/(sqrt(float(PI)+log(eij)*log(eij))); //generalized damping - Cummins 2010
                      //const float gn=0.08f*sqrt(nu_mass*sqrt(float(Dp)/2)/((taup1+taup2)/2)); //generalized damping - Lemieux 2008
                float rep=kn*pow(over_lap,1.5f);
                      float fn=rep-gn*pow(over_lap,0.25f)*vn;                   
                      acep1.x+=(fn*nx); acep1.y+=(fn*ny); acep1.z+=(fn*nz); //-Force is applied in the normal between the particles
                      //tangencial
                      float dvxt=dvx-vn*nx, dvyt=dvy-vn*ny, dvzt=dvz-vn*nz; //Vji_t
                      float vt=sqrt(dvxt*dvxt + dvyt*dvyt + dvzt*dvzt);
                      float tx=0, ty=0, tz=0; //Tang vel unit vector
                if(vt!=0){ tx=dvxt/vt; ty=dvyt/vt; tz=dvzt/vt; }
                      float ft_elast=2*(kn*float(DemDtForce)-gn)*vt/7;   //Elastic frictional string -->  ft_elast=2*(kn*fdispl-gn*vt)/7; fdispl=dtforce*vt;
                const float kfric_ij=(kfricp1+kfricp2)/2;
                      float ft=kfric_ij*fn*tanh(8*vt);  //Coulomb
                      ft=(ft<ft_elast? ft: ft_elast);   //not above yield criteria, visco-elastic model
                      acep1.x+=(ft*tx); acep1.y+=(ft*ty); acep1.z+=(ft*tz);
                    }
            }
          }
        }
      }
      //-Almacena resultados.
      if(acep1.x||acep1.y||acep1.z){
        ace[p1]=ace[p1]+acep1;
        const int th=omp_get_thread_num();
        if(demdtth[th*STRIDE_OMP]<demdtp1)demdtth[th*STRIDE_OMP]=demdtp1;
      }
    }
  }
  //-Actualiza viscdt con el valor maximo de viscdt y demdt*.
  float demdt=demdtth[0];
  for(int th=1;th<OmpThreads;th++)if(demdt<demdtth[th*STRIDE_OMP])demdt=demdtth[th*STRIDE_OMP];
  if(viscdt<demdt)viscdt=demdt;
}


//==============================================================================
/// Computes sub-particle stress tensor (Tau) for SPS turbulence model.   
//==============================================================================
void JSphCpu::ComputeSpsTau(unsigned n,unsigned pini,const tfloat4 *velrhop,const tsymatrix3f *gradvel,tsymatrix3f *tau)const{
  const int pfin=int(pini+n);
  #ifdef _WITHOMP
    #pragma omp parallel for schedule (static)
  #endif
  for(int p=int(pini);p<pfin;p++){
    const tsymatrix3f gradvel=SpsGradvelc[p];
    const float pow1=gradvel.xx*gradvel.xx + gradvel.yy*gradvel.yy + gradvel.zz*gradvel.zz;
    const float prr=pow1+pow1 + gradvel.xy*gradvel.xy + gradvel.xz*gradvel.xz + gradvel.yz*gradvel.yz;
    const float visc_sps=SpsSmag*sqrt(prr);
    const float div_u=gradvel.xx+gradvel.yy+gradvel.zz;
    const float sps_k=(2.0f/3.0f)*visc_sps*div_u;
    const float sps_blin=SpsBlin*prr;
    const float sumsps=-(sps_k+sps_blin);
    const float twovisc_sps=(visc_sps+visc_sps);
    const float one_rho2=1.0f/velrhop[p].w;   
    tau[p].xx=one_rho2*(twovisc_sps*gradvel.xx +sumsps);
    tau[p].xy=one_rho2*(visc_sps   *gradvel.xy);
    tau[p].xz=one_rho2*(visc_sps   *gradvel.xz);
    tau[p].yy=one_rho2*(twovisc_sps*gradvel.yy +sumsps);
    tau[p].yz=one_rho2*(visc_sps   *gradvel.yz);
    tau[p].zz=one_rho2*(twovisc_sps*gradvel.zz +sumsps);
  }
}


//==============================================================================
// Seleccion de parametros template para Interaction_ForcesFluidT.
//==============================================================================
template<bool psimple,TpFtMode ftmode,bool lamsps,TpDeltaSph tdelta,bool shift> void JSphCpu::Interaction_ForcesT
  (unsigned np,unsigned npb,unsigned npbok
  ,tuint3 ncells,const unsigned *begincell,tuint3 cellmin,const unsigned *dcell
  ,const tdouble3 *pos,const tfloat3 *pspos,const tfloat4 *velrhop,const word *code,const unsigned *idp
  ,const float *press
  ,float &viscdt,float* ar,tfloat3 *ace,float *delta
  ,tsymatrix3f *spstau,tsymatrix3f *spsgradvel
  ,TpShifting tshifting,tfloat3 *shiftpos,float *shiftdetect)
{
  const unsigned npf=np-npb;
  const tint4 nc=TInt4(int(ncells.x),int(ncells.y),int(ncells.z),int(ncells.x*ncells.y));
  const tint3 cellzero=TInt3(cellmin.x,cellmin.y,cellmin.z);
  const unsigned cellfluid=nc.w*nc.z+1;
  const int hdiv=(CellMode==CELLMODE_H? 2: 1);
  
  if(npf){
    //-Interaccion Fluid-Fluid
    InteractionForcesFluid<psimple,ftmode,lamsps,tdelta,shift>               (npf,npb,nc,hdiv,cellfluid,Visco                 ,begincell,cellzero,dcell,spstau,spsgradvel,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,tshifting,shiftpos,shiftdetect);
    //-Interaccion Fluid-Bound
    InteractionForcesFluid<psimple,ftmode,lamsps,tdelta,shift> (npf,npb,nc,hdiv,0        ,Visco*ViscoBoundFactor,begincell,cellzero,dcell,spstau,spsgradvel,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,tshifting,shiftpos,shiftdetect);

    //-Interaccion DEM Floating-Bound & Floating-Floating //(DEM)
    if(USE_DEM)InteractionForcesDEM<psimple> (CaseNfloat,nc,hdiv,cellfluid,begincell,cellzero,dcell,FtRidp,DemObjs,pos,pspos,velrhop,code,idp,viscdt,ace);

    //-Interaccion con viscosidad Laminar+SPS
    if(lamsps){
      //InteractionForcesSps<psimple,ftmode> (npf,npb,nc,hdiv,cellfluid,Visco,begincell,cellzero,dcell,pos,pspos,velrhop,code,idp,spstau,spsgradvel,ace);
      ComputeSpsTau(npf,npb,velrhop,spsgradvel,spstau);
    }
  }
  if(npbok){
    //-Interaccion Bound-Fluid
    InteractionForcesBound      <psimple,ftmode> (npbok,0,nc,hdiv,cellfluid,begincell,cellzero,dcell,pos,pspos,velrhop,code,idp,viscdt,ar);
  }
  //-Para simulaciones 2D anula siempre la 2� componente
  if(Simulate2D && npf)for(unsigned p=Npb;p<Np;p++)ace[p].y=0;
}

//==============================================================================
// Seleccion de parametros template para Interaction_ForcesX.
//==============================================================================
void JSphCpu::Interaction_Forces(unsigned np,unsigned npb,unsigned npbok
  ,tuint3 ncells,const unsigned *begincell,tuint3 cellmin,const unsigned *dcell
  ,const tdouble3 *pos,const tfloat4 *velrhop,const unsigned *idp,const word *code
  ,const float *press
  ,float &viscdt,float* ar,tfloat3 *ace,float *delta
  ,tsymatrix3f *spstau,tsymatrix3f *spsgradvel
  ,tfloat3 *shiftpos,float *shiftdetect)
{
  tfloat3 *pspos=NULL;
  const bool psimple=false;
  if(!WithFloating){                   const TpFtMode ftmode=FTMODE_None;
    if(TShifting){                     const bool tshift=true;
      if(TVisco==VISCO_LaminarSPS){    const bool lamsps=true;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }else{                           const bool lamsps=false;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }
    }else{                             const bool tshift=false;
      if(TVisco==VISCO_LaminarSPS){    const bool lamsps=true;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }else{                           const bool lamsps=false;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }
    }
  }else if(!UseDEM){                   const TpFtMode ftmode=FTMODE_Sph;
    if(TShifting){                     const bool tshift=true;
      if(TVisco==VISCO_LaminarSPS){    const bool lamsps=true;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }else{                           const bool lamsps=false;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }
    }else{                             const bool tshift=false;
      if(TVisco==VISCO_LaminarSPS){    const bool lamsps=true;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }else{                           const bool lamsps=false;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }
    }
  }else{                               const TpFtMode ftmode=FTMODE_Dem;
    if(TShifting){                     const bool tshift=true;
      if(TVisco==VISCO_LaminarSPS){    const bool lamsps=true;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }else{                           const bool lamsps=false;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }
    }else{                             const bool tshift=false;
      if(TVisco==VISCO_LaminarSPS){    const bool lamsps=true;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }else{                           const bool lamsps=false;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }
    }
  }
}

//==============================================================================
// Seleccion de parametros template para Interaction_ForcesX.
//==============================================================================
void JSphCpu::InteractionSimple_Forces(unsigned np,unsigned npb,unsigned npbok
  ,tuint3 ncells,const unsigned *begincell,tuint3 cellmin,const unsigned *dcell
  ,const tfloat3 *pspos,const tfloat4 *velrhop,const unsigned *idp,const word *code
  ,const float *press
  ,float &viscdt,float* ar,tfloat3 *ace,float *delta
  ,tsymatrix3f *spstau,tsymatrix3f *spsgradvel
  ,tfloat3 *shiftpos,float *shiftdetect)
{
  tdouble3 *pos=NULL;
  const bool psimple=true;
  if(!WithFloating){                   const TpFtMode ftmode=FTMODE_None;
    if(TShifting){                     const bool tshift=true;
      if(TVisco==VISCO_LaminarSPS){    const bool lamsps=true;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }else{                           const bool lamsps=false;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }
    }else{                             const bool tshift=false;
      if(TVisco==VISCO_LaminarSPS){    const bool lamsps=true;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }else{                           const bool lamsps=false;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }
    }
  }else if(!UseDEM){                   const TpFtMode ftmode=FTMODE_Sph;
    if(TShifting){                     const bool tshift=true;
      if(TVisco==VISCO_LaminarSPS){    const bool lamsps=true;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }else{                           const bool lamsps=false;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }
    }else{                             const bool tshift=false;
      if(TVisco==VISCO_LaminarSPS){    const bool lamsps=true;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }else{                           const bool lamsps=false;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }
    }
  }else{                               const TpFtMode ftmode=FTMODE_Dem;
    if(TShifting){                     const bool tshift=true;
      if(TVisco==VISCO_LaminarSPS){    const bool lamsps=true;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }else{                           const bool lamsps=false;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }
    }else{                             const bool tshift=false;
      if(TVisco==VISCO_LaminarSPS){    const bool lamsps=true;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }else{                           const bool lamsps=false;
        if(TDeltaSph==DELTA_None)      Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_None,tshift>       (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_Dynamic)   Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_Dynamic,tshift>    (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
        if(TDeltaSph==DELTA_DynamicExt)Interaction_ForcesT<psimple,ftmode,lamsps,DELTA_DynamicExt,tshift> (np,npb,npbok,ncells,begincell,cellmin,dcell,pos,pspos,velrhop,code,idp,press,viscdt,ar,ace,delta,spstau,spsgradvel,TShifting,shiftpos,shiftdetect);
      }
    }
  }
}

//==============================================================================
// Actualiza pos, dcell y code a partir del desplazamiento indicado.
// El valor de outrhop indica si esta fuera de los limites de densidad.
// Comprueba los limites en funcion de MapRealPosMin y MapRealSize esto es valido
// para single-cpu pq DomRealPos y MapRealPos son iguales. Para multi-cpu seria 
// necesario marcar las particulas q salgan del dominio sin salir del mapa.
//==============================================================================
void JSphCpu::UpdatePos(tdouble3 rpos,double movx,double movy,double movz
  ,bool outrhop,unsigned p,tdouble3 *pos,unsigned *cell,word *code)const
{
  //-Comprueba validez del desplazamiento.
  bool outmove=(fabs(float(movx))>MovLimit || fabs(float(movy))>MovLimit || fabs(float(movz))>MovLimit);
  //-Aplica desplazamiento.
  rpos.x+=movx; rpos.y+=movy; rpos.z+=movz;
  //-Comprueba limites del dominio reales.
  double dx=rpos.x-MapRealPosMin.x;
  double dy=rpos.y-MapRealPosMin.y;
  double dz=rpos.z-MapRealPosMin.z;
  bool out=(dx!=dx || dy!=dy || dz!=dz || dx<0 || dy<0 || dz<0 || dx>=MapRealSize.x || dy>=MapRealSize.y || dz>=MapRealSize.z);
  //-Ajusta posicion segun condiciones periodicas y vuelve a comprobar los limites del dominio.
  if(PeriActive && out){
    bool xperi=((PeriActive&1)!=0),yperi=((PeriActive&2)!=0),zperi=((PeriActive&4)!=0);
    if(xperi){
      if(dx<0)             { dx-=PeriXinc.x; dy-=PeriXinc.y; dz-=PeriXinc.z; }
      if(dx>=MapRealSize.x){ dx+=PeriXinc.x; dy+=PeriXinc.y; dz+=PeriXinc.z; }
    }
    if(yperi){
      if(dy<0)             { dx-=PeriYinc.x; dy-=PeriYinc.y; dz-=PeriYinc.z; }
      if(dy>=MapRealSize.y){ dx+=PeriYinc.x; dy+=PeriYinc.y; dz+=PeriYinc.z; }
    }
    if(zperi){
      if(dz<0)             { dx-=PeriZinc.x; dy-=PeriZinc.y; dz-=PeriZinc.z; }
      if(dz>=MapRealSize.z){ dx+=PeriZinc.x; dy+=PeriZinc.y; dz+=PeriZinc.z; }
    }
    bool outx=!xperi && (dx<0 || dx>=MapRealSize.x);
    bool outy=!yperi && (dy<0 || dy>=MapRealSize.y);
    bool outz=!zperi && (dz<0 || dz>=MapRealSize.z);
    out=(outx||outy||outz);
    rpos=TDouble3(dx,dy,dz)+MapRealPosMin;
  }
  //-Guarda posicion actualizada.
  pos[p]=rpos;
  //-Guarda celda y check.
  if(outrhop || outmove || out){//-Particle out
    word rcode=code[p];
    if(outrhop)rcode=CODE_SetOutRhop(rcode);
    else if(out)rcode=CODE_SetOutPos(rcode);
    else rcode=CODE_SetOutMove(rcode);
    code[p]=rcode;
    cell[p]=0xFFFFFFFF;
  }
  else{//-Particle in
    if(PeriActive){
      dx=rpos.x-DomPosMin.x;
      dy=rpos.y-DomPosMin.y;
      dz=rpos.z-DomPosMin.z;
    }
    unsigned cx=unsigned(dx/Scell),cy=unsigned(dy/Scell),cz=unsigned(dz/Scell);
    cell[p]=PC__Cell(DomCellCode,cx,cy,cz);
  }
}

//==============================================================================
// Calcula nuevos valores de posicion, velocidad y densidad para el fluido (usando Verlet).
//==============================================================================
template<bool shift> void JSphCpu::ComputeVerletVarsFluid(const tfloat4 *velrhop1,const tfloat4 *velrhop2,double dt,double dt2
  ,tdouble3 *pos,unsigned *dcell,word *code,tfloat4 *velrhopnew)const
{
  const double dt205=0.5*dt*dt;
  const int pini=int(Npb),pfin=int(Np),npf=int(Np-Npb);
  #ifdef _WITHOMP
    #pragma omp parallel for schedule (static) if(npf>LIMIT_COMPUTESTEP_OMP)
  #endif
  for(int p=pini;p<pfin;p++){
    //-Calcula densidad.
    const float rhopnew=float(double(velrhop2[p].w)+dt2*Arc[p]);
    if(!WithFloating || CODE_GetType(code[p])==CODE_TYPE_FLUID){//-Particulas: Fluid
      //-Calcula desplazamiento y actualiza posicion.
      double dx=double(velrhop1[p].x)*dt + double(Acec[p].x)*dt205;
      double dy=double(velrhop1[p].y)*dt + double(Acec[p].y)*dt205;
      double dz=double(velrhop1[p].z)*dt + double(Acec[p].z)*dt205;
      if(shift){
        //if(Idpc[p]==5759)Log->Printf("%u_SHIFT_COMPUTE> p[%u] px:%f pz:%f",Nstep,p,ShiftPosc[p].x,ShiftPosc[p].z); 
        dx+=double(ShiftPosc[p].x);
        dy+=double(ShiftPosc[p].y);
        dz+=double(ShiftPosc[p].z);
      }
      bool outrhop=(rhopnew<RhopOutMin||rhopnew>RhopOutMax);
      UpdatePos(pos[p],dx,dy,dz,outrhop,p,pos,dcell,code);
      //-Actualiza velocidad y densidad.
      velrhopnew[p].x=float(double(velrhop2[p].x)+double(Acec[p].x)*dt2);
      velrhopnew[p].y=float(double(velrhop2[p].y)+double(Acec[p].y)*dt2);
      velrhopnew[p].z=float(double(velrhop2[p].z)+double(Acec[p].z)*dt2);
      velrhopnew[p].w=rhopnew;
    }
    else{//-Particulas: Floating
      velrhopnew[p]=velrhop1[p];
      velrhopnew[p].w=(rhopnew<RHOPCERO? RHOPCERO: rhopnew); //-Evita q las floating absorvan a las fluidas.
    }
  }
}

//==============================================================================
// Calcula nuevos valores de densidad y pone velocidad a cero para el contorno 
// (fixed+moving, no floating).
//==============================================================================
void JSphCpu::ComputeVelrhopBound(const tfloat4* velrhopold,double armul,tfloat4* velrhopnew)const{
  const int npb=int(Npb);
  #ifdef _WITHOMP
    #pragma omp parallel for schedule (static) if(npb>LIMIT_COMPUTESTEP_OMP)
  #endif
  for(int p=0;p<npb;p++){
    const float rhopnew=float(double(velrhopold[p].w)+armul*Arc[p]);
    velrhopnew[p]=TFloat4(0,0,0,(rhopnew<RHOPCERO? RHOPCERO: rhopnew));//-Evita q las boundary absorvan a las fluidas.
  }
}

//==============================================================================
// Actualizacion de particulas segun fuerzas y dt usando Verlet.
//==============================================================================
void JSphCpu::ComputeVerlet(double dt){
  TmcStart(Timers,TMC_SuComputeStep);
  VerletStep++;
  if(VerletStep<VerletSteps){
    const double twodt=dt+dt;
    if(TShifting)ComputeVerletVarsFluid<true>  (Velrhopc,VelrhopM1c,dt,twodt,Posc,Dcellc,Codec,VelrhopM1c);
    else         ComputeVerletVarsFluid<false> (Velrhopc,VelrhopM1c,dt,twodt,Posc,Dcellc,Codec,VelrhopM1c);
    ComputeVelrhopBound(VelrhopM1c,twodt,VelrhopM1c);
  }
  else{
    if(TShifting)ComputeVerletVarsFluid<true>  (Velrhopc,Velrhopc,dt,dt,Posc,Dcellc,Codec,VelrhopM1c);
    else         ComputeVerletVarsFluid<false> (Velrhopc,Velrhopc,dt,dt,Posc,Dcellc,Codec,VelrhopM1c);
    ComputeVelrhopBound(Velrhopc,dt,VelrhopM1c);
    VerletStep=0;
  }
  //-Los nuevos valores se calculan en VelrhopM1c.
  swap(Velrhopc,VelrhopM1c);     //-Intercambia Velrhopc y VelrhopM1c.
  TmcStop(Timers,TMC_SuComputeStep);
}

//==============================================================================
// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Predictor.
//==============================================================================
void JSphCpu::ComputeSymplecticPre(double dt){
  if(TShifting)ComputeSymplecticPreT<true> (dt);
  else         ComputeSymplecticPreT<false>(dt);
}
//==============================================================================
template<bool shift> void JSphCpu::ComputeSymplecticPreT(double dt){
  TmcStart(Timers,TMC_SuComputeStep);
  //-Asigna memoria a variables Pre.
  PosPrec=ArraysCpu->ReserveDouble3();
  VelrhopPrec=ArraysCpu->ReserveFloat4();
  //-Cambia datos a variables Pre para calcular nuevos datos.
  swap(PosPrec,Posc);         //Es decir... PosPre[] <= Pos[]
  swap(VelrhopPrec,Velrhopc); //Es decir... VelrhopPre[] <= Velrhop[]
  //-Calcula nuevos datos de particulas.
  const double dt05=dt*.5;
  
  //-Calcula nueva densidad para el contorno y copia velocidad.
  const int npb=int(Npb);
  #ifdef _WITHOMP
    #pragma omp parallel for schedule (static) if(npb>LIMIT_COMPUTESTEP_OMP)
  #endif
  for(int p=0;p<npb;p++){
    const tfloat4 vr=VelrhopPrec[p];
    const float rhopnew=float(double(vr.w)+dt05*Arc[p]);
    Velrhopc[p]=TFloat4(vr.x,vr.y,vr.z,(rhopnew<RHOPCERO? RHOPCERO: rhopnew));//-Evita q las boundary absorvan a las fluidas.
  }

  //-Calcula nuevos datos del fluido.
  const int np=int(Np);
  #ifdef _WITHOMP
    #pragma omp parallel for schedule (static) if(np>LIMIT_COMPUTESTEP_OMP)
  #endif
  for(int p=npb;p<np;p++){
    //-Calcula densidad.
    const float rhopnew=float(double(VelrhopPrec[p].w)+dt05*Arc[p]);
    if(!WithFloating || CODE_GetType(Codec[p])==CODE_TYPE_FLUID){//-Particulas: Fluid
      //-Calcula desplazamiento y actualiza posicion.
      double dx=double(VelrhopPrec[p].x)*dt05;
      double dy=double(VelrhopPrec[p].y)*dt05;
      double dz=double(VelrhopPrec[p].z)*dt05;
      if(shift){
        dx+=double(ShiftPosc[p].x);
        dy+=double(ShiftPosc[p].y);
        dz+=double(ShiftPosc[p].z);
      }
      bool outrhop=(rhopnew<RhopOutMin||rhopnew>RhopOutMax);
      UpdatePos(PosPrec[p],dx,dy,dz,outrhop,p,Posc,Dcellc,Codec);
      //-Actualiza velocidad y densidad.
      Velrhopc[p].x=float(double(VelrhopPrec[p].x)+double(Acec[p].x)* dt05);
      Velrhopc[p].y=float(double(VelrhopPrec[p].y)+double(Acec[p].y)* dt05);
      Velrhopc[p].z=float(double(VelrhopPrec[p].z)+double(Acec[p].z)* dt05);
      Velrhopc[p].w=rhopnew;
    }
    else{//-Particulas: Floating
      Velrhopc[p]=VelrhopPrec[p];
      Velrhopc[p].w=(rhopnew<RHOPCERO? RHOPCERO: rhopnew); //-Evita q las floating absorvan a las fluidas.
      //-Copia posicion.
      Posc[p]=PosPrec[p];
      //Dcellc[p]=0xFFFFFFFF;  //-Para que de error sino se actualiza en RunFloating(). (No es compatible con FtPause).
    }
  }

  //-Copia posicion anterior del contorno.
  memcpy(Posc,PosPrec,sizeof(tdouble3)*Npb);

  TmcStop(Timers,TMC_SuComputeStep);
}

//==============================================================================
// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Corrector.
//==============================================================================
void JSphCpu::ComputeSymplecticCorr(double dt){
  if(TShifting)ComputeSymplecticCorrT<true> (dt);
  else         ComputeSymplecticCorrT<false>(dt);
}
//==============================================================================
template<bool shift> void JSphCpu::ComputeSymplecticCorrT(double dt){
  TmcStart(Timers,TMC_SuComputeStep);
  
  //-Calcula rhop de contorno y vel igual a cero.
  const int npb=int(Npb);
  #ifdef _WITHOMP
    #pragma omp parallel for schedule (static) if(npb>LIMIT_COMPUTESTEP_OMP)
  #endif
  for(int p=0;p<npb;p++){
    const double epsilon_rdot=(-double(Arc[p])/double(Velrhopc[p].w))*dt;
    const float rhopnew=float(double(VelrhopPrec[p].w) * (2.-epsilon_rdot)/(2.+epsilon_rdot));
    Velrhopc[p]=TFloat4(0,0,0,(rhopnew<RHOPCERO? RHOPCERO: rhopnew));//-Evita q las boundary absorvan a las fluidas.
  }

  //-Calcula datos de fluido.
  const double dt05=dt*.5;
  const int np=int(Np);
  #ifdef _WITHOMP
    #pragma omp parallel for schedule (static) if(np>LIMIT_COMPUTESTEP_OMP)
  #endif
  for(int p=npb;p<np;p++){
    const double epsilon_rdot=(-double(Arc[p])/double(Velrhopc[p].w))*dt;
    const float rhopnew=float(double(VelrhopPrec[p].w) * (2.-epsilon_rdot)/(2.+epsilon_rdot));
    if(!WithFloating || CODE_GetType(Codec[p])==CODE_TYPE_FLUID){//-Particulas: Fluid
      //-Actualiza velocidad y densidad.
      Velrhopc[p].x=float(double(VelrhopPrec[p].x) + double(Acec[p].x) * dt); 
      Velrhopc[p].y=float(double(VelrhopPrec[p].y) + double(Acec[p].y) * dt); 
      Velrhopc[p].z=float(double(VelrhopPrec[p].z) + double(Acec[p].z) * dt); 
      Velrhopc[p].w=rhopnew;
      //-Calcula desplazamiento y actualiza posicion.
      double dx=(double(VelrhopPrec[p].x)+double(Velrhopc[p].x)) * dt05; 
      double dy=(double(VelrhopPrec[p].y)+double(Velrhopc[p].y)) * dt05; 
      double dz=(double(VelrhopPrec[p].z)+double(Velrhopc[p].z)) * dt05;
      if(shift){
        dx+=double(ShiftPosc[p].x);
        dy+=double(ShiftPosc[p].y);
        dz+=double(ShiftPosc[p].z);
      }
      bool outrhop=(rhopnew<RhopOutMin||rhopnew>RhopOutMax);
      UpdatePos(PosPrec[p],dx,dy,dz,outrhop,p,Posc,Dcellc,Codec);
    }
    else{//-Particulas: Floating
      Velrhopc[p]=VelrhopPrec[p];
      Velrhopc[p].w=(rhopnew<RHOPCERO? RHOPCERO: rhopnew); //-Evita q las floating absorvan a las fluidas.
      //-Copia posicion.
      Posc[p]=PosPrec[p];
      //Dcellc[p]=0xFFFFFFFF;  //-Para que de error sino se actualiza en RunFloating(). (No es compatible con FtPause).
    }
  }

  //-Libera memoria asignada a variables Pre en ComputeSymplecticPre().
  ArraysCpu->Free(PosPrec);      PosPrec=NULL;
  ArraysCpu->Free(VelrhopPrec);  VelrhopPrec=NULL;
  TmcStop(Timers,TMC_SuComputeStep);
}

//==============================================================================
// Calcula un Dt variable.
//==============================================================================
double JSphCpu::DtVariable(bool final){
  //-dt1 depends on force per unit mass.
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
void JSphCpu::RunShifting(double dt){
  const double coeftfs=(Simulate2D? 2.0: 3.0)-ShiftTFS;
  const int pini=int(Npb),pfin=int(Np),npf=int(Np-Npb);
  #ifdef _WITHOMP
    #pragma omp parallel for schedule (static) if(npf>LIMIT_COMPUTELIGHT_OMP)
  #endif
  for(int p=pini;p<pfin;p++){
    double vx=double(Velrhopc[p].x);
    double vy=double(Velrhopc[p].y);
    double vz=double(Velrhopc[p].z);
    double umagn=double(ShiftCoef)*double(H)*sqrt(vx*vx+vy*vy+vz*vz)*dt;
    if(ShiftDetectc){
      if(ShiftDetectc[p]<ShiftTFS)umagn=0;
      else umagn*=(double(ShiftDetectc[p])-ShiftTFS)/coeftfs;
    }
    if(ShiftPosc[p].x==FLT_MAX)umagn=0; //-Anula shifting por proximidad del contorno.
    ShiftPosc[p]=ToTFloat3(ToTDouble3(ShiftPosc[p])*umagn);
  }
  //if(PartNstep+1==Nstep)SaveVtkGradDiv(Part-1);
}

//==============================================================================
// Calcula posicion de particulas segun idp[]. Cuando no la encuentra es UINT_MAX.
// Cuando periactive es False sumpone que no hay particulas duplicadas (periodicas)
// y todas son CODE_NORMAL.
//==============================================================================
void JSphCpu::CalcRidp(bool periactive,unsigned np,unsigned pini,unsigned idini,unsigned idfin,const word *code,const unsigned *idp,unsigned *ridp)const{
  //-Asigna valores UINT_MAX
  const unsigned nsel=idfin-idini;
  memset(ridp,255,sizeof(unsigned)*nsel); 

  //-Calcula posicion segun id.
  const unsigned pfin=pini+np;
  if(periactive){//-Calcula posicion segun id comprobando que las particulas son normales (no periodicas).
    for(unsigned p=pini;p<pfin;p++){
      const unsigned id=idp[p];
      if(idini<=id && id<idfin){
        if(CODE_GetSpecialValue(code[p])==CODE_NORMAL)ridp[id-idini]=p;
      }
    }
  }
  else{//-Calcula posicion segun id suponiendo que todas las particulas son normales (no periodicas).
    for(unsigned p=pini;p<pfin;p++){
      const unsigned id=idp[p];
      if(idini<=id && id<idfin)ridp[id-idini]=p;
    }
  }
}

//==============================================================================
// Aplica un movimiento lineal a un conjunto de particulas.
//==============================================================================
void JSphCpu::MoveLinBound(unsigned np,unsigned ini,const tdouble3 &mvpos,const tfloat3 &mvvel
  ,const unsigned *ridp,tdouble3 *pos,unsigned *dcell,tfloat4 *velrhop,word *code)const
{
  const unsigned fin=ini+np;
  for(unsigned id=ini;id<fin;id++){
    const unsigned pid=RidpMove[id];  //printf("id:%d -> pid:%d\n",id,pid);
    //Log->Printf("id:%d -> pid:%d",id,pid);
    if(pid!=UINT_MAX){
      //Log->Printf("  mvpos:(%f,%f,%f)",mvpos.x,mvpos.y,mvpos.z);
      //Log->Printf("  mvvel:(%f,%f,%f)",mvvel.x,mvvel.y,mvvel.z);
      UpdatePos(pos[pid],mvpos.x,mvpos.y,mvpos.z,false,pid,pos,dcell,code);
      velrhop[pid].x=mvvel.x;  velrhop[pid].y=mvvel.y;  velrhop[pid].z=mvvel.z;
    }
  }
}

//==============================================================================
// Aplica un movimiento matricial a un conjunto de particulas.
//==============================================================================
void JSphCpu::MoveMatBound(unsigned np,unsigned ini,tmatrix4d m,double dt
  ,const unsigned *ridpmv,tdouble3 *pos,unsigned *dcell,tfloat4 *velrhop,word *code)const
{
  const unsigned fin=ini+np;
  for(unsigned id=ini;id<fin;id++){
    const unsigned pid=RidpMove[id];
    if(pid!=UINT_MAX){
      tdouble3 ps=pos[pid];
      tdouble3 ps2=MatrixMulPoint(m,ps);
      if(Simulate2D)ps2.y=ps.y;
      const double dx=ps2.x-ps.x, dy=ps2.y-ps.y, dz=ps2.z-ps.z;
      UpdatePos(ps,dx,dy,dz,false,pid,pos,dcell,code);
      velrhop[pid].x=float(dx/dt);  velrhop[pid].y=float(dy/dt);  velrhop[pid].z=float(dz/dt);
    }
  }
}

//==============================================================================
// Procesa movimiento de boundary particles
//==============================================================================
void JSphCpu::RunMotion(double stepdt){
  const char met[]="RunMotion";
  TmcStart(Timers,TMC_SuMotion);

  unsigned nmove=0;
  if(Motion->ProcesTime(TimeStep+MotionTimeMod,stepdt)){
    nmove=Motion->GetMovCount();
    //{ char cad[256]; sprintf(cad,"----RunMotion[%u]>  nmove:%u",Nstep,nmove); Log->Print(cad); }
    if(nmove){
      CalcRidp(PeriActive!=0,Npb,0,CaseNfixed,CaseNfixed+CaseNmoving,Codec,Idpc,RidpMove);
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
          MoveLinBound(np,pini,mvsimple,mvvel,RidpMove,Posc,Dcellc,Velrhopc,Codec);
        }
        else{//-Movimiento con matriz
          const unsigned pini=MotionObjBegin[ref]-CaseNfixed,np=MotionObjBegin[ref+1]-MotionObjBegin[ref];
          mvmatrix=OrderCode(mvmatrix);
          MoveMatBound(np,pini,mvmatrix,stepdt,RidpMove,Posc,Dcellc,Velrhopc,Codec);
        }
      }
      BoundChanged=true;
    }
  }
  //-Procesa otros modos de motion.
  if(WaveGen){
    if(!nmove)CalcRidp(PeriActive!=0,Npb,0,CaseNfixed,CaseNfixed+CaseNmoving,Codec,Idpc,RidpMove);
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
        //Log->Printf("------> mvsimple:(%f,%f,%f) stepdt:%f",mvsimple.x,mvsimple.y,mvsimple.z,stepdt);
        MoveLinBound(nparts,idbegin-CaseNfixed,mvsimple,mvvel,RidpMove,Posc,Dcellc,Velrhopc,Codec);
      }
      else{
        mvmatrix=OrderCode(mvmatrix);
        MoveMatBound(nparts,idbegin-CaseNfixed,mvmatrix,stepdt,RidpMove,Posc,Dcellc,Velrhopc,Codec);
      }
    }
  }
  TmcStop(Timers,TMC_SuMotion);
}

//==============================================================================
// Ajusta variables de particulas floating body
//==============================================================================
void JSphCpu::InitFloating(){
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
}

//==============================================================================
// Muestra los temporizadores activos.
//==============================================================================
void JSphCpu::ShowTimers(bool onlyfile){
  JLog2::TpMode_Out mode=(onlyfile? JLog2::Out_File: JLog2::Out_ScrFile);
  Log->Print("\n[CPU Timers]",mode);
  if(!SvTimers)Log->Print("none",mode);
  else for(unsigned c=0;c<TimerGetCount();c++)if(TimerIsActive(c))Log->Print(TimerToText(c),mode);
}

//==============================================================================
// Devuelve string con nombres y valores de los timers activos.
//==============================================================================
void JSphCpu::GetTimersInfo(std::string &hinfo,std::string &dinfo)const{
  for(unsigned c=0;c<TimerGetCount();c++)if(TimerIsActive(c)){
    hinfo=hinfo+";"+TimerGetName(c);
    dinfo=dinfo+";"+fun::FloatStr(TimerGetValue(c)/1000.f);
  }
}


