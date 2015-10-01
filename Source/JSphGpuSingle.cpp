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

#include "JSphGpuSingle.h"
#include "JCellDivGpuSingle.h"
#include "JArraysGpu.h"
#include "Functions.h"
#include "JXml.h"
#include "JSphMotion.h"
#include "JPartsLoad4.h"
#include "JSphVisco.h"
#include "JWaveGen.h"
#include "JSphGpu_ker.h"
#include "JPtxasInfo.h"

#include "JFormatFiles2.h"

using namespace std;
//==============================================================================
// Constructor.
//==============================================================================
JSphGpuSingle::JSphGpuSingle():JSphGpu(false){
  ClassName="JSphGpuSingle";
  CellDivSingle=NULL;
  PartsLoaded=NULL;
}

//==============================================================================
// Destructor.
//==============================================================================
JSphGpuSingle::~JSphGpuSingle(){
  delete CellDivSingle; CellDivSingle=NULL;
  delete PartsLoaded;   PartsLoaded=NULL;
}

//==============================================================================
// Devuelve la memoria reservada en cpu.
//==============================================================================
llong JSphGpuSingle::GetAllocMemoryCpu()const{  
  llong s=JSphGpu::GetAllocMemoryCpu();
  //Reservada en otros objetos
  if(CellDivSingle)s+=CellDivSingle->GetAllocMemoryCpu();
  if(PartsLoaded)s+=PartsLoaded->GetAllocMemory();
  return(s);
}

//==============================================================================
// Devuelve la memoria reservada en gpu.
//==============================================================================
llong JSphGpuSingle::GetAllocMemoryGpu()const{  
  llong s=JSphGpu::GetAllocMemoryGpu();
  //Reservada en otros objetos
  if(CellDivSingle)s+=CellDivSingle->GetAllocMemoryGpu();
  return(s);
}

//==============================================================================
// Devuelve la memoria gpu reservada o usada para particulas.
//==============================================================================
llong JSphGpuSingle::GetMemoryGpuNp()const{
  llong s=JSphGpu::GetAllocMemoryGpu();
  //Reservada en otros objetos
  if(CellDivSingle)s+=CellDivSingle->GetAllocMemoryGpuNp();
  return(s);
}

//==============================================================================
// Devuelve la memoria gpu reservada o usada para celdas.
//==============================================================================
llong JSphGpuSingle::GetMemoryGpuNct()const{
  llong s=CellDivSingle->GetAllocMemoryGpuNct();
  return(CellDivSingle->GetAllocMemoryGpuNct());
}

//==============================================================================
// Actualiza los valores maximos de memory, particles y cells.
//==============================================================================
void JSphGpuSingle::UpdateMaxValues(){
  MaxParticles=max(MaxParticles,Np);
  if(CellDivSingle)MaxCells=max(MaxCells,CellDivSingle->GetNct());
  llong m=GetAllocMemoryCpu();
  MaxMemoryCpu=max(MaxMemoryCpu,m);
  m=GetAllocMemoryGpu();
  MaxMemoryGpu=max(MaxMemoryGpu,m);
}

//==============================================================================
// Carga la configuracion de ejecucion.
//==============================================================================
void JSphGpuSingle::LoadConfig(JCfgRun *cfg){
  const char met[]="LoadConfig";
  PtxasFile=cfg->PtxasFile;
  //-Carga configuracion basica general
  JSph::LoadConfig(cfg);
  //-Checks compatibility of selected options.
  if(RenCorrection && UseDEM)RunException(met,"Ren correction is not implemented with Floatings-DEM.");
  Log->Print("**Special case configuration is loaded");
}

//==============================================================================
// Carga particulas del caso a procesar.
//==============================================================================
void JSphGpuSingle::LoadCaseParticles(){
  Log->Print("Loading initial state of particles...");
  PartsLoaded=new JPartsLoad4;
  PartsLoaded->LoadParticles(DirCase,CaseName,PartBegin,PartBeginDir);
  PartsLoaded->CheckConfig(CaseNp,CaseNfixed,CaseNmoving,CaseNfloat,CaseNfluid,PeriX,PeriY,PeriZ);
  Log->Printf("Loaded particles: %u",PartsLoaded->GetCount());
  //-Recupera informacion de las particulas cargadas.
  Simulate2D=PartsLoaded->GetSimulate2D();
  if(Simulate2D&&PeriY)RunException("LoadCaseParticles","Can not use periodic conditions in Y with 2D simulations");
  CasePosMin=PartsLoaded->GetCasePosMin();
  CasePosMax=PartsLoaded->GetCasePosMax();

  //-Calcula limites reales de la simulacion.
  if(PartsLoaded->MapSizeLoaded())PartsLoaded->GetMapSize(MapRealPosMin,MapRealPosMax);
  else{
    PartsLoaded->CalculeLimits(double(H)*BORDER_MAP,Dp/2.,PeriX,PeriY,PeriZ,MapRealPosMin,MapRealPosMax);
    ResizeMapLimits();
  }
  if(PartBegin){
    PartBeginTimeStep=PartsLoaded->GetPartBeginTimeStep();
    PartBeginTotalNp=PartsLoaded->GetPartBeginTotalNp();
  }
  Log->Print(string("MapRealPos(final)=")+fun::Double3gRangeStr(MapRealPosMin,MapRealPosMax));
  MapRealSize=MapRealPosMax-MapRealPosMin;
  Log->Print("**Initial state of particles is loaded");

  //-Configura limites de ejes periodicos
  if(PeriX)PeriXinc.x=-MapRealSize.x;
  if(PeriY)PeriYinc.y=-MapRealSize.y;
  if(PeriZ)PeriZinc.z=-MapRealSize.z;
  //-Calcula limites de simulacion con bordes periodicos.
  Map_PosMin=MapRealPosMin; Map_PosMax=MapRealPosMax;
  float dosh=float(H*2);
  if(PeriX){ Map_PosMin.x=Map_PosMin.x-dosh;  Map_PosMax.x=Map_PosMax.x+dosh; }
  if(PeriY){ Map_PosMin.y=Map_PosMin.y-dosh;  Map_PosMax.y=Map_PosMax.y+dosh; }
  if(PeriZ){ Map_PosMin.z=Map_PosMin.z-dosh;  Map_PosMax.z=Map_PosMax.z+dosh; }
  Map_Size=Map_PosMax-Map_PosMin;
}

//==============================================================================
// Configuracion del dominio actual.
//==============================================================================
void JSphGpuSingle::ConfigDomain(){
  const char* met="ConfigDomain";
  //-Calcula numero de particulas.
  Np=PartsLoaded->GetCount(); Npb=CaseNpb; NpbOk=Npb;
  //-Reserva memoria fija para moving y floating.
  AllocGpuMemoryFixed();
  //-Reserva memoria en Gpu para particulas.
  AllocGpuMemoryParticles(Np,0);
  //-Reserva memoria en Cpu.
  AllocCpuMemoryParticles(Np);

  //-Copia datos de particulas.
  memcpy(AuxPos,PartsLoaded->GetPos(),sizeof(tdouble3)*Np);
  memcpy(Idp,PartsLoaded->GetIdp(),sizeof(unsigned)*Np);
  memcpy(Velrhop,PartsLoaded->GetVelRhop(),sizeof(tfloat4)*Np);

  //-Calcula radio de floatings.
  if(CaseNfloat && PeriActive!=0)CalcFloatingRadius(Np,AuxPos,Idp);

  //-Carga code de particulas.
  LoadCodeParticles(Np,Idp,Code);

  //-Libera memoria de PartsLoaded.
  delete PartsLoaded; PartsLoaded=NULL;
  //-Aplica configuracion de CellOrder.
  ConfigCellOrder(CellOrder,Np,AuxPos,Velrhop);

  //-Configura division celdas.
  ConfigCellDivision();
  //-Establece dominio de simulacion local dentro de Map_Cells y calcula DomCellCode.
  SelecDomain(TUint3(0,0,0),Map_Cells);
  //-Calcula celda inicial de particulas y comprueba si hay excluidas inesperadas.
  LoadDcellParticles(Np,Code,AuxPos,Dcell);

  //-Sube datos de particulas a la GPU.
  ReserveBasicArraysGpu();
  for(unsigned p=0;p<Np;p++){ Posxy[p]=TDouble2(AuxPos[p].x,AuxPos[p].y); Posz[p]=AuxPos[p].z; }
  ParticlesDataUp(Np);
  //-Sube constantes a la GPU.
  ConstantDataUp();

  //-Crea objeto para divide en Gpu y selecciona un cellmode valido.
  CellDivSingle=new JCellDivGpuSingle(Stable,FtCount!=0,PeriActive,CellOrder,CellMode,Scell,Map_PosMin,Map_PosMax,Map_Cells,CaseNbound,CaseNfixed,CaseNpb,Log,DirOut);
  CellDivSingle->DefineDomain(DomCellCode,DomCelIni,DomCelFin,DomPosMin,DomPosMax);
  ConfigCellDiv((JCellDivGpu*)CellDivSingle);

  ConfigBlockSizes(false,PeriActive!=0);
  ConfigSaveData(0,1,"");

  //-Reordena particulas por celda.
  BoundChanged=true;
  RunCellDivide(true);
}

//==============================================================================
// Redimensiona el espacio reservado para particulas en CPU y GPU midiendo el
// tiempo consumido con TMG_SuResizeNp. Al terminar actualiza el divide.
//==============================================================================
void JSphGpuSingle::ResizeParticlesSize(unsigned newsize,float oversize,bool updatedivide){
  TmgStart(Timers,TMG_SuResizeNp);
  newsize+=(oversize>0? unsigned(oversize*newsize): 0);
  FreeCpuMemoryParticles();
  CellDivSingle->FreeMemoryGpu();
  ResizeGpuMemoryParticles(newsize);
  AllocCpuMemoryParticles(newsize);
  TmgStop(Timers,TMG_SuResizeNp);
  if(updatedivide)RunCellDivide(true);
}

//==============================================================================
// Crea particulas duplicadas de condiciones periodicas.
// Crea nuevas particulas periodicas y marca las viejas para ignorarlas.
// Las nuevas periodicas se situan a partir del Np de entrada, primero las NpbPer
// de contorno y despues las NpfPer fluidas. El Np de salida contiene tambien las
// nuevas periodicas.
//==============================================================================
void JSphGpuSingle::RunPeriodic(){
  const char met[]="RunPeriodic";
  TmgStart(Timers,TMG_SuPeriodic);
  //-Guarda numero de periodicas actuales.
  NpfPerM1=NpfPer;
  NpbPerM1=NpbPer;
  //-Marca periodicas actuales para ignorar.
  cusph::PeriodicIgnore(Np,Codeg);
  //-Crea las nuevas periodicas.
  const unsigned npb0=Npb;
  const unsigned npf0=Np-Npb;
  const unsigned np0=Np;
  NpbPer=NpfPer=0;
  BoundChanged=true;
  for(unsigned ctype=0;ctype<2;ctype++){//-0:bound, 1:fluid+floating.
    //-Calcula rango de particulas a examinar (bound o fluid).
    const unsigned pini=(ctype? npb0: 0);
    const unsigned num= (ctype? npf0: npb0);
    //-Busca periodicas en cada eje (X, Y e Z).
    for(unsigned cper=0;cper<3;cper++)if((cper==0 && PeriActive&1) || (cper==1 && PeriActive&2) || (cper==2 && PeriActive&4)){
      tdouble3 perinc=(cper==0? PeriXinc: (cper==1? PeriYinc: PeriZinc));
      //-Primero busca en la lista de periodicas nuevas y despues en la lista inicial de particulas (necesario para periodicas en mas de un eje).
      for(unsigned cblock=0;cblock<2;cblock++){//-0:periodicas nuevas, 1:particulas originales
        const unsigned nper=(ctype? NpfPer: NpbPer); //-Numero de periodicas nuevas del tipo a procesar.
        const unsigned pini2=(cblock? pini: Np-nper);
        const unsigned num2= (cblock? num:  nper);
        //-Repite la busqueda si la memoria disponible resulto insuficiente y hubo que aumentarla.
        bool run=true;
        while(run && num2){
          //-Reserva memoria para crear lista de particulas periodicas.
          unsigned* listpg=ArraysGpu->ReserveUint();
          unsigned nmax=GpuParticlesSize-1; //-Numero maximo de particulas que caben en la lista.
          //-Genera lista de nuevas periodicas.
          if(Np>=0x80000000)RunException(met,"The number of particles is too big.");//-Pq el ultimo bit se usa para marcar el sentido en que se crea la nueva periodica.
          unsigned count=cusph::PeriodicMakeList(num2,pini2,Stable,nmax,Map_PosMin,Map_PosMax,perinc,Posxyg,Poszg,Codeg,listpg);
          //-Redimensiona memoria para particulas si no hay espacio suficiente y repite el proceso de busqueda.
          if(count>nmax || count+Np>GpuParticlesSize){
            ArraysGpu->Free(listpg); listpg=NULL;
            TmgStop(Timers,TMG_SuPeriodic);
            ResizeParticlesSize(Np+count,PERIODIC_OVERMEMORYNP,false);
            TmgStart(Timers,TMG_SuPeriodic);
          }
          else{
            run=false;
            //-Crea nuevas particulas periodicas duplicando las particulas de la lista.
            if(TStep==STEP_Verlet)cusph::PeriodicDuplicateVerlet(count,Np,DomCells,perinc,listpg,Idpg,Codeg,Dcellg,Posxyg,Poszg,Velrhopg,SpsTaug,VelrhopM1g);
            if(TStep==STEP_Symplectic){
              if((PosxyPreg || PoszPreg || VelrhopPreg) && (!PosxyPreg || !PoszPreg || !VelrhopPreg))RunException(met,"Symplectic data is invalid.") ;
              cusph::PeriodicDuplicateSymplectic(count,Np,DomCells,perinc,listpg,Idpg,Codeg,Dcellg,Posxyg,Poszg,Velrhopg,SpsTaug,PosxyPreg,PoszPreg,VelrhopPreg);
            }
            //-Libera lista y actualiza numero de particulas.
            ArraysGpu->Free(listpg); listpg=NULL;
            Np+=count;
            //-Actualiza numero de periodicas nuevas.
            if(!ctype)NpbPer+=count;
            else NpfPer+=count;
          }
        }
      }
    }
  }
  TmgStop(Timers,TMG_SuPeriodic);
  CheckCudaError(met,"Failed in creation of periodic particles.");
}

//==============================================================================
// Ejecuta divide de particulas en celdas.
//==============================================================================
void JSphGpuSingle::RunCellDivide(bool updateperiodic){
  const char met[]="RunCellDivide";
  //-Crea nuevas particulas periodicas y marca las viejas para ignorarlas.
  if(updateperiodic && PeriActive)RunPeriodic();

  //-Inicia Divide.
  CellDivSingle->Divide(Npb,Np-Npb-NpbPer-NpfPer,NpbPer,NpfPer,BoundChanged,Dcellg,Codeg,Timers,Posxyg,Poszg,Idpg);

  //-Ordena datos de particulas
  TmgStart(Timers,TMG_NlSortData);
  {
    unsigned* idpg=ArraysGpu->ReserveUint();
    word*     codeg=ArraysGpu->ReserveWord();
    unsigned* dcellg=ArraysGpu->ReserveUint();
    double2*  posxyg=ArraysGpu->ReserveDouble2();
    double*   poszg=ArraysGpu->ReserveDouble();
    float4*   velrhopg=ArraysGpu->ReserveFloat4();
    CellDivSingle->SortBasicArrays(Idpg,Codeg,Dcellg,Posxyg,Poszg,Velrhopg,idpg,codeg,dcellg,posxyg,poszg,velrhopg);
    swap(Idpg,idpg);           ArraysGpu->Free(idpg);
    swap(Codeg,codeg);         ArraysGpu->Free(codeg);
    swap(Dcellg,dcellg);       ArraysGpu->Free(dcellg);
    swap(Posxyg,posxyg);       ArraysGpu->Free(posxyg);
    swap(Poszg,poszg);         ArraysGpu->Free(poszg);
    swap(Velrhopg,velrhopg);   ArraysGpu->Free(velrhopg);
  }
  if(TStep==STEP_Verlet){
    float4* velrhopg=ArraysGpu->ReserveFloat4();
    CellDivSingle->SortDataArrays(VelrhopM1g,velrhopg);
    swap(VelrhopM1g,velrhopg);   ArraysGpu->Free(velrhopg);
  }
  else if(TStep==STEP_Symplectic && (PosxyPreg || PoszPreg || VelrhopPreg)){//En realidad solo es necesario en el divide del corrector, no en el predictor???
    if(!PosxyPreg || !PoszPreg || !VelrhopPreg)RunException(met,"Symplectic data is invalid.") ;
    double2* posxyg=ArraysGpu->ReserveDouble2();
    double* poszg=ArraysGpu->ReserveDouble();
    float4* velrhopg=ArraysGpu->ReserveFloat4();
    CellDivSingle->SortDataArrays(PosxyPreg,PoszPreg,VelrhopPreg,posxyg,poszg,velrhopg);
    swap(PosxyPreg,posxyg);      ArraysGpu->Free(posxyg);
    swap(PoszPreg,poszg);        ArraysGpu->Free(poszg);
    swap(VelrhopPreg,velrhopg);  ArraysGpu->Free(velrhopg);
  }
  if(TVisco==VISCO_LaminarSPS){
    tsymatrix3f *spstaug=ArraysGpu->ReserveSymatrix3f();
    CellDivSingle->SortDataArrays(SpsTaug,spstaug);
    swap(SpsTaug,spstaug);  ArraysGpu->Free(spstaug);
  }

  //-Recupera datos del divide.
  Np=CellDivSingle->GetNpFinal();
  Npb=CellDivSingle->GetNpbFinal();
  NpbOk=Npb-CellDivSingle->GetNpbIgnore();
  //-Recupera posiciones de floatings.
  if(CaseNfloat)cusph::CalcRidp(PeriActive!=0,Np-Npb,Npb,CaseNpb,CaseNpb+CaseNfloat,Codeg,Idpg,FtRidpg);
  TmgStop(Timers,TMG_NlSortData);

  //-Gestion de particulas excluidas (contorno y fluid).
  TmgStart(Timers,TMG_NlOutCheck);
  unsigned npout=CellDivSingle->GetNpOut();
  if(npout){
    ParticlesDataDown(npout,Np,true,true,false);
    CellDivSingle->CheckParticlesOut(npout,Idp,AuxPos,AuxRhop,Code);
    AddParticlesOut(npout,Idp,AuxPos,AuxVel,AuxRhop,CellDivSingle->GetNpfOutRhop(),CellDivSingle->GetNpfOutMove());
  }
  TmgStop(Timers,TMG_NlOutCheck);
  BoundChanged=false;
}

//==============================================================================
// Aplica correccion de Ren a la presion y densidad del contorno.
//==============================================================================
void JSphGpuSingle::RunRenCorrection(){
  //-Calcula presion en contorno a partir de fluido.
  float *presskf=ArraysGpu->ReserveFloat();
  cusph::Interaction_Ren(Psimple,WithFloating,CellMode,NpbOk,CellDivSingle->GetNcells(),CellDivSingle->GetBeginCell(),CellDivSingle->GetCellDomainMin(),Dcellg,Posxyg,Poszg,PsPospressg,Velrhopg,Codeg,Idpg,FtoMasspg,Gravity,presskf);
  //-Recalcula valores de presion y densidad en contorno segun RenBeta.
  cusph::ComputeRenPress(Psimple,NpbOk,RenCorrection,presskf,Velrhopg,PsPospressg);
  ArraysGpu->Free(presskf); presskf=NULL;
}

//==============================================================================
// Interaccion para el calculo de fuerzas.
//==============================================================================
void JSphGpuSingle::Interaction_Forces(TpInter tinter){
  const char met[]="Interaction_Forces";
  PreInteraction_Forces(tinter);
  TmgStart(Timers,TMG_CfForces);
  if(RenCorrection)RunRenCorrection();

  const bool lamsps=(TVisco==VISCO_LaminarSPS);
  const unsigned bsfluid=BlockSizes.forcesfluid;
  const unsigned bsbound=BlockSizes.forcesbound;

  //-Interaccion Fluid-Fluid/Bound & Bound-Fluid
  cusph::Interaction_Forces(Psimple,WithFloating,UseDEM,lamsps,TDeltaSph,CellMode,Visco*ViscoBoundFactor,Visco,bsbound,bsfluid,Np,Npb,NpbOk,CellDivSingle->GetNcells(),CellDivSingle->GetBeginCell(),CellDivSingle->GetCellDomainMin(),Dcellg,Posxyg,Poszg,PsPospressg,Velrhopg,Codeg,Idpg,FtoMasspg,SpsTaug,SpsGradvelg,ViscDtg,Arg,Aceg,Deltag,TShifting,ShiftPosg,ShiftDetectg,Simulate2D);
  
  //if(1){  //dbg
  //  unsigned *idph=new unsigned[Np];
  //  tfloat4 *pospressh=new tfloat4[Np];
  //  tfloat4 *velrhoph=new tfloat4[Np];
  //  float *arh=new float[Np];
  //  float3 *aceh=new float3[Np];
  //  //float3 *shifth=new float3[Np];
  //  cudaMemcpy(idph,Idpg,sizeof(unsigned)*Np,cudaMemcpyDeviceToHost);
  //  cudaMemcpy(pospressh,PsPospressg,sizeof(float4)*Np,cudaMemcpyDeviceToHost);
  //  cudaMemcpy(velrhoph,Velrhopg,sizeof(float4)*Np,cudaMemcpyDeviceToHost);
  //  cudaMemcpy(arh,Arg,sizeof(float)*Np,cudaMemcpyDeviceToHost);
  //  cudaMemcpy(aceh,Aceg,sizeof(float3)*Np,cudaMemcpyDeviceToHost);
  //  //cudaMemcpy(shifth,ShiftPosg,sizeof(float3)*Np,cudaMemcpyDeviceToHost);
  //  unsigned idsel=508;
  //  for(unsigned p=0;p<Np;p++)if(idph[p]==idsel){
  //    Log->Printf("%u> particle[%u]> idp:%u  ar:%f  ace:(%f,%f,%f) ",Nstep,p,idph[p],arh[p],aceh[p].x,aceh[p].y,aceh[p].z);
  //  }
  //  Log->Print(" ");
  //  delete[] idph;
  //  delete[] pospressh;
  //  delete[] velrhoph;
  //  delete[] arh;
  //  delete[] aceh;
  //}

  //-Interaccion DEM Floating-Bound & Floating-Floating //(DEM)
  if(UseDEM)cusph::Interaction_ForcesDem(Psimple,CellMode,BlockSizes.forcesdem,CaseNfloat,CellDivSingle->GetNcells(),CellDivSingle->GetBeginCell(),CellDivSingle->GetCellDomainMin(),Dcellg,FtRidpg,DemDatag,float(DemDtForce),Posxyg,Poszg,PsPospressg,Velrhopg,Codeg,Idpg,ViscDtg,Aceg);

  //-Calculo de Tau para Laminar+SPS
  if(lamsps)cusph::ComputeSpsTau(Np,Npb,SpsSmag,SpsBlin,Velrhopg,SpsGradvelg,SpsTaug);

  if(Deltag)cusph::AddDelta(Np-Npb,Deltag+Npb,Arg+Npb);//-Añade correccion de Delta-SPH a Arg[].
  CheckCudaError(met,"Failed while executing kernels of interaction.");

  //-Calculates maximum value of ViscDt.
  if(Np)ViscDtMax=cusph::ReduMaxFloat(Np,0,ViscDtg,CellDivSingle->GetAuxMem(cusph::ReduMaxFloatSize(Np)));
  //-Calculates maximum value of Ace.
  AceMax=ComputeAceMax(ViscDtg); 

  TmgStop(Timers,TMG_CfForces);
  CheckCudaError(met,"Failed in reduction of viscdt.");
}

//==============================================================================
// Devuelve valor maximo de (ace.x^2 + ace.y^2 + ace.z^2) a partir de Acec[].
//==============================================================================
double JSphGpuSingle::ComputeAceMax(float *auxmem){
  float acemax=0;
  if(!PeriActive)cusph::ComputeAceMod(Np-Npb,Aceg+Npb,auxmem);//-Sin condiciones periodicas.
  else cusph::ComputeAceMod(Np-Npb,Codeg+Npb,Aceg+Npb,auxmem);//-Con condiciones periodicas ignora las particulas periodicas.
  if(Np-Npb)acemax=cusph::ReduMaxFloat(Np-Npb,0,auxmem,CellDivSingle->GetAuxMem(cusph::ReduMaxFloatSize(Np-Npb)));
  return(sqrt(double(acemax)));
}

//==============================================================================
// Genera fichero VTK con datos de particulas (debug)
//==============================================================================
void JSphGpuSingle::DgSaveVtk(std::string fname,unsigned num,bool svace){
  const unsigned np=Np;
  //-Reserva memoria.
  tfloat3 *pos=new tfloat3[np];
  tfloat3 *ace=NULL;
  float *ar=NULL;
  float *viscdt=NULL;
  tfloat3 *velcor=NULL;
  if(svace){
    ace=new tfloat3[np];
    ar=new float[np];
    viscdt=new float[np];
  }
  //-Recupera datos de particulas.
  ParticlesDataDown(np,0,true,true,false);
  for(unsigned p=0;p<np;p++)pos[p]=ToTFloat3(AuxPos[p]);
  if(ace)cudaMemcpy(ace,Aceg,sizeof(tfloat3)*np,cudaMemcpyDeviceToHost);
  if(ar)cudaMemcpy(ar,Arg,sizeof(float)*np,cudaMemcpyDeviceToHost);
  if(viscdt)cudaMemcpy(viscdt,ViscDtg,sizeof(float)*np,cudaMemcpyDeviceToHost);
  //-Define campos de VTK.
  JFormatFiles2::StScalarData fields[10];
  unsigned nfields=0;
  if(Idp){     fields[nfields]=JFormatFiles2::DefineField("Id"  ,JFormatFiles2::UInt32  ,1,Idp);     nfields++; }
  if(Code){    fields[nfields]=JFormatFiles2::DefineField("Code",JFormatFiles2::UShort16,1,Code);    nfields++; }
  if(AuxVel){  fields[nfields]=JFormatFiles2::DefineField("Vel" ,JFormatFiles2::Float32 ,3,AuxVel);  nfields++; }
  if(AuxRhop){ fields[nfields]=JFormatFiles2::DefineField("Rhop",JFormatFiles2::Float32 ,1,AuxRhop); nfields++; }
  if(ace){     fields[nfields]=JFormatFiles2::DefineField("Ace" ,JFormatFiles2::Float32 ,3,ace);     nfields++; }
  if(ar){      fields[nfields]=JFormatFiles2::DefineField("Ar"  ,JFormatFiles2::Float32 ,1,ar);      nfields++; }
  if(viscdt){  fields[nfields]=JFormatFiles2::DefineField("Visc",JFormatFiles2::Float32 ,1,viscdt);  nfields++; }
  //-Genera fichero.
  JFormatFiles2::SaveVtk(DirOut+fun::FileNameSec(fname+".vtk",num),np,pos,nfields,fields);
  JFormatFiles2::SaveCsv(DirOut+fun::FileNameSec(fname+".csv",num),np,pos,nfields,fields);
  //-Libera memoria.
  delete[] pos;    pos=NULL;
  delete[] ace;    ace=NULL;
  delete[] ar;     ar=NULL;
  delete[] viscdt; viscdt=NULL;
}

//==============================================================================
// Realiza interaccion y actualizacion de particulas segun las fuerzas 
// calculadas en la interaccion usando Verlet.
//==============================================================================
double JSphGpuSingle::ComputeStep_Ver(){
  //DgSaveVtk("_A_PreInter",Nstep,false);
  Interaction_Forces(INTER_Forces);     //-Interaccion
  const double dt=DtVariable(true);     //-Calcula nuevo dt
  DemDtForce=dt;                        //(DEM)
  if(TShifting)RunShifting(dt);         //-Shifting
  ComputeVerlet(dt);                    //-Actualiza particulas usando Verlet
  if(CaseNfloat)RunFloating(dt,false); //-Gestion de floating bodies
  PosInteraction_Forces();              //-Libera memoria de interaccion
  return(dt);
}

//==============================================================================
// Realiza interaccion y actualizacion de particulas segun las fuerzas 
// calculadas en la interaccion usando Symplectic.
//==============================================================================
double JSphGpuSingle::ComputeStep_Sym(){
  const double dt=DtPre;
  //-Predictor
  //-----------
  DemDtForce=dt*0.5f;                     //(DEM)
  Interaction_Forces(INTER_Forces);       //-Interaccion
  const double ddt_p=DtVariable(false);   //-Calcula dt del predictor
  if(TShifting)RunShifting(dt*.5);        //-Shifting
  ComputeSymplecticPre(dt);               //-Aplica Symplectic-Predictor a las particulas
  if(CaseNfloat)RunFloating(dt*.5,true); //-Gestion de floating bodies
  PosInteraction_Forces();                //-Libera memoria de interaccion
  //-Corrector
  //-----------
  DemDtForce=dt;                          //(DEM)
  RunCellDivide(true);
  Interaction_Forces(INTER_ForcesCorr);   //-Interaccion
  const double ddt_c=DtVariable(true);    //-Calcula dt del corrector
  if(TShifting)RunShifting(dt);           //-Shifting
  ComputeSymplecticCorr(dt);              //-Aplica Symplectic-Corrector a las particulas
  if(CaseNfloat)RunFloating(dt,false);   //-Gestion de floating bodies
  PosInteraction_Forces();                //-Libera memoria de interaccion

  DtPre=min(ddt_p,ddt_c);                 //-Calcula el dt para el siguiente ComputeStep
  return(dt);
}

//==============================================================================
// Procesa floating objects.
//==============================================================================
void JSphGpuSingle::RunFloating(double dt,bool predictor){
  const char met[]="RunFloating";
  //RunFloating_Old(dt,predictor);
  if(TimeStep>=FtPause){//-Se usa >= pq si FtPause es cero en symplectic-predictor no entraria.
    TmgStart(Timers,TMG_SuFloating);
    //-Gets positions of floating particles.
    //cusph::CalcRidp(PeriActive!=0,Np-Npb,Npb,CaseNpb,CaseNpb+CaseNfloat,Codeg,Idpg,FtRidpg); Se hace en RunCellDivide
    //-Calcula fuerzas sobre floatings.
    cusph::FtCalcForces(PeriActive!=0,FtCount,Gravity,FtoDatag,FtoMasspg,FtoCenterg,FtRidpg,Posxyg,Poszg,Aceg,FtoForcesg);
    //-Aplica movimiento sobre floatings.
    cusph::FtUpdate(PeriActive!=0,predictor,Simulate2D,FtCount,dt,Gravity,FtoDatag,FtRidpg,FtoForcesg,FtoCenterg,FtoVelg,FtoOmegag,Posxyg,Poszg,Dcellg,Velrhopg,Codeg);
    TmgStop(Timers,TMG_SuFloating);
  }
}

//==============================================================================
// Procesa floating objects (codigo viejo).
// Atencion:Para usar esta version de RunFloating se debe comentar la recuperacion
// de datos de floatings en Save().
//==============================================================================
void JSphGpuSingle::RunFloating_Old(double dt2,bool predictor){
  if(TimeStep>=FtPause){//-Se usa >= pq si FtPause es cero en symplectic-predictor no entraria.
    TmgStart(Timers,TMG_SuFloating);
    //cusph::CalcRidp(PeriActive!=0,Np-Npb,Npb,CaseNpb,CaseNpb+CaseNfloat,Codeg,Idpg,FtRidpg); Se hace en RunCellDivide
    for(unsigned cf=0;cf<FtCount;cf++){
      StFloatingData *fobj=FtObjs+cf;
      const float fradius=fobj->radius;
      tdouble3 fcenter=fobj->center;
      //-Computes traslational and rotational velocities.
      tfloat3 face,fomegavel;
      tmatrix3f inert=TMatrix3f(0,0,0,0,0,0,0,0,0);
      {
        float3 *resultg=(float3 *)CellDivSingle->GetAuxMem(15);
        cusph::FtCalcOmega(PeriActive!=0,fobj->count,fobj->begin-CaseNpb,Gravity,fobj->massp,fcenter,fobj->radius,FtRidpg,Posxyg,Poszg,Aceg,resultg);
        tfloat3 result[5];
        cudaMemcpy(result,resultg,sizeof(float3)*5,cudaMemcpyDeviceToHost);
        face=result[0];
        fomegavel=result[1];
        inert=TMatrix3f(result[2].x,result[2].y,result[2].z,result[3].x,result[3].y,result[3].z,result[4].x,result[4].y,result[4].z);
      }
      //Log->Printf("%u__FT>> face:%s fomegavel:%s",Nstep,fun::Float3Str(face,"%f,%f,%f").c_str(),fun::Float3Str(fomegavel,"%f,%f,%f").c_str());

      face.x=(face.x+fobj->mass*Gravity.x)/fobj->mass;
      face.y=(face.y+fobj->mass*Gravity.y)/fobj->mass;
      face.z=(face.z+fobj->mass*Gravity.z)/fobj->mass;
      //Log->Printf("%u__FT_1>> face:%s ",Nstep,fun::Float3Str(face,"%f,%f,%f").c_str());

      //-Recomputes values of floating.
      tfloat3 fvel=fobj->fvel;
      //Log->Printf("%u__FT_1>> center:%s ",Nstep,fun::Double3Str(center,"%f,%f,%.10f").c_str());
      //Log->Printf("%u__FT_1>> fvel:%s ",Nstep,fun::Float3Str(fvel,"%f,%f,%.10f").c_str());

      //calculates the inverse of the intertia matrix to compute the I^-1 * L= W
      tmatrix3f invinert=TMatrix3f(0,0,0,0,0,0,0,0,0);
      float detiner=(inert.a11*inert.a22*inert.a33+inert.a12*inert.a23*inert.a31+inert.a21*inert.a32*inert.a13-(inert.a31*inert.a22*inert.a13+inert.a21*inert.a12*inert.a33+inert.a23*inert.a32*inert.a11));
      if(detiner){
        invinert.a11=(inert.a22*inert.a33-inert.a23*inert.a32)/detiner;
        invinert.a12=-(inert.a12*inert.a33-inert.a13*inert.a32)/detiner;
        invinert.a13=(inert.a12*inert.a23-inert.a13*inert.a22)/detiner;
        invinert.a21=-(inert.a21*inert.a33-inert.a23*inert.a31)/detiner;
        invinert.a22=(inert.a11*inert.a33-inert.a13*inert.a31)/detiner;
        invinert.a23=-(inert.a11*inert.a23-inert.a13*inert.a21)/detiner;
        invinert.a31=(inert.a21*inert.a32-inert.a22*inert.a31)/detiner;
        invinert.a32=-(inert.a11*inert.a32-inert.a12*inert.a31)/detiner;
        invinert.a33=(inert.a11*inert.a22-inert.a12*inert.a21)/detiner;
      }

      tfloat3 omega=TFloat3(0);
      //correct omega formulation
      omega.x=fomegavel.x*invinert.a11+fomegavel.y*invinert.a12+fomegavel.z*invinert.a13;
      omega.y=fomegavel.x*invinert.a21+fomegavel.y*invinert.a22+fomegavel.z*invinert.a23;
      omega.z=fomegavel.x*invinert.a31+fomegavel.y*invinert.a32+fomegavel.z*invinert.a33;
      //Log->Printf("%u__FT_1>> detiner:%f omega:%s",Nstep,detiner,fun::Float3Str(omega,"%f,%f,%f").c_str());
      tfloat3 fomega;
      fomega.x=float(dt2*omega.x+fobj->fomega.x);
      fomega.y=float(dt2*omega.y+fobj->fomega.y);
      fomega.z=float(dt2*omega.z+fobj->fomega.z);
      //Log->Printf("%u__FT_1>> fomega:%s ",Nstep,fun::Float3Str(fomega,"%f,%f,%f").c_str());
      if(Simulate2D){ face.y=0; fomega.x=0; fomega.z=0; fvel.y=0; }
      fcenter.x+=dt2*fvel.x;
      fcenter.y+=dt2*fvel.y;
      fcenter.z+=dt2*fvel.z;
      //Log->Printf("%u__FT_2>> center:%s ",Nstep,fun::Double3Str(center,"%f,%f,%f").c_str());
      fvel.x=float(dt2*face.x+fvel.x);
      fvel.y=float(dt2*face.y+fvel.y);
      fvel.z=float(dt2*face.z+fvel.z);
      //Log->Printf("%u__FT_2>> fvel:%s ",Nstep,fun::Float3Str(fvel,"%f,%f,%f").c_str());

      //-Updates floating particles.
      cusph::FtUpdate(PeriActive!=0,predictor,fobj->count,fobj->begin-CaseNpb,dt2,fcenter,fradius,fvel,fomega,FtRidpg,Posxyg,Poszg,Dcellg,Velrhopg,Codeg);

      //-Stores data.
      if(!predictor){
        fobj->center=(PeriActive? UpdatePeriodicPos(fcenter): fcenter);
        fobj->fvel=fvel;
        fobj->fomega=fomega;
      }
    }
    TmgStop(Timers,TMG_SuFloating);
  }
}

//==============================================================================
// Inicia proceso de simulacion.
//==============================================================================
void JSphGpuSingle::Run(std::string appname,JCfgRun *cfg,JLog2 *log){
  const char* met="Run";
  if(!cfg||!log)return;
  AppName=appname; Log=log;

  //-Seleccion de GPU
  //-------------------
  SelecDevice(cfg->GpuId);

  //-Configura timers
  //-------------------
  TmgCreation(Timers,cfg->SvTimers);
  TmgStart(Timers,TMG_Init);
  if(cfg->SvTimersStep>0){
    TimersStep=new JTimersStep(cfg->DirOut,cfg->SvTimersStep,0,0);
    for(unsigned ct=0;ct<TimerGetCount();ct++)if(TimerIsActive(ct))TimersStep->AddTimer(TimerGetName(ct),TimerGetPtrValue(ct));
  }

  //-Carga de parametros y datos de entrada
  //-----------------------------------------
  LoadConfig(cfg);
  LoadCaseParticles();
  ConfigConstants(Simulate2D);
  ConfigDomain();
  ConfigRunMode("Single-Gpu");

  //-Inicializacion de variables de ejecucion
  //-------------------------------------------
  InitRun();
  UpdateMaxValues();
  PrintAllocMemory(GetAllocMemoryCpu(),GetAllocMemoryGpu());
  SaveData(); 
  TmgResetValues(Timers);
  TmgStop(Timers,TMG_Init);
  PartNstep=-1; Part++;

  //-Bucle principal
  //------------------
  bool partoutstop=false;
  TimerSim.Start();
  TimerPart.Start();
  Log->Print(string("\n[Initialising simulation (")+RunCode+")  "+fun::GetDateTime()+"]");
  PrintHeadPart();
  while(TimeStep<TimeMax){
    if(ViscoTime)Visco=ViscoTime->GetVisco(float(TimeStep));
    double stepdt=ComputeStep();
    if(PartDtMin>stepdt)PartDtMin=stepdt; if(PartDtMax<stepdt)PartDtMax=stepdt;
    if(CaseNmoving)RunMotion(stepdt);
    RunCellDivide(true);
    TimeStep+=stepdt;
    partoutstop=(Np<NpMinimum || !Np);
    if((TimeStep-TimeStepIni)-TimePart*((Part-PartIni)-1)>=TimePart || partoutstop){
      if(partoutstop){
        Log->Print("\n**** Particles OUT limit reached...\n");
        TimeMax=TimeStep;
      }
      SaveData();
      Part++;
      PartNstep=Nstep;
      TimeStepM1=TimeStep;
      TimerPart.Start();
    }
    UpdateMaxValues();
    Nstep++;
    if(TimersStep&&TimersStep->Check(float(TimeStep)))SaveTimersStep(Np,Npb,NpbOk,CellDivSingle->GetNct());
    //if(Nstep>=1)break;
  }
  TimerSim.Stop(); TimerTot.Stop();

  //-Fin de simulacion
  //--------------------
  FinishRun(partoutstop);
}

//==============================================================================
// Genera los ficheros de salida de datos
//==============================================================================
void JSphGpuSingle::SaveData(){
  const bool save=(SvData!=SDAT_None&&SvData!=SDAT_Info);
  const unsigned npsave=Np-NpbPer-NpfPer; //-Resta las periodicas si las hubiera.
  //-Recupera datos de particulas en GPU.
  if(save){
    TmgStart(Timers,TMG_SuDownData);
    unsigned npnormal=ParticlesDataDown(Np,0,false,true,PeriActive!=0);
    if(npnormal!=npsave)RunException("SaveData","The number of particles is invalid.");
    TmgStop(Timers,TMG_SuDownData);
  }
  //-Recupera datos de floatings en GPU.
  if(FtCount){
    TmgStart(Timers,TMG_SuDownData);
    cudaMemcpy(FtoCenter,FtoCenterg,sizeof(double3)*FtCount,cudaMemcpyDeviceToHost);
    for(unsigned cf=0;cf<FtCount;cf++)FtObjs[cf].center=FtoCenter[cf];
    tfloat3 *aux=(tfloat3 *)FtoCenter;
    cudaMemcpy(aux,FtoVelg,sizeof(float3)*FtCount,cudaMemcpyDeviceToHost);
    for(unsigned cf=0;cf<FtCount;cf++)FtObjs[cf].fvel=aux[cf];
    cudaMemcpy(aux,FtoOmegag,sizeof(float3)*FtCount,cudaMemcpyDeviceToHost);
    for(unsigned cf=0;cf<FtCount;cf++)FtObjs[cf].fomega=aux[cf];
    TmgStop(Timers,TMG_SuDownData);
  }
  //-Reune informacion adicional.
  TmgStart(Timers,TMG_SuSavePart);
  StInfoPartPlus infoplus;
  memset(&infoplus,0,sizeof(StInfoPartPlus));
  if(SvData&SDAT_Info){
    infoplus.nct=CellDivSingle->GetNct();
    infoplus.npbin=NpbOk;
    infoplus.npbout=Npb-NpbOk;
    infoplus.npf=Np-Npb;
    infoplus.npbper=NpbPer;
    infoplus.npfper=NpfPer;
    infoplus.memorycpualloc=this->GetAllocMemoryCpu();
    infoplus.gpudata=true;
    infoplus.memorynctalloc=infoplus.memorynctused=GetMemoryGpuNct();
    infoplus.memorynpalloc=infoplus.memorynpused=GetMemoryGpuNp();
    TimerSim.Stop();
    infoplus.timesim=TimerSim.GetElapsedTimeD()/1000.;
  }
  //-Graba datos de particulas.
  const tdouble3 vdom[2]={OrderDecode(CellDivSingle->GetDomainLimits(true)),OrderDecode(CellDivSingle->GetDomainLimits(false))};
  JSph::SaveData(npsave,Idp,AuxPos,AuxVel,AuxRhop,1,vdom,&infoplus);
  //-Graba informacion de ejecucion.
  if(TimersStep)TimersStep->SaveData();
  TmgStop(Timers,TMG_SuSavePart);
}

//==============================================================================
// Muestra y graba resumen final de ejecucion.
//==============================================================================
void JSphGpuSingle::FinishRun(bool stop){
  float tsim=TimerSim.GetElapsedTimeF()/1000.f,ttot=TimerTot.GetElapsedTimeF()/1000.f;
  if(TimersStep)TimersStep->SaveData();
  JSph::ShowResume(stop,tsim,ttot,true,"");
  string hinfo=";RunMode",dinfo=string(";")+RunMode;
  if(SvTimers){
    ShowTimers();
    GetTimersInfo(hinfo,dinfo);
  }
  Log->Print(" ");
  if(SvRes)SaveRes(tsim,ttot,hinfo,dinfo);
}




//==============================================================================
// Graba fichero vtk con datos de las particulas.
//==============================================================================
void DgSaveVtkParticlesGpuX(std::string filename,int numfile,unsigned pini,unsigned pfin,unsigned cellcode
  ,const double2 *posxyg,const double *poszg,const unsigned *idpg,const unsigned *dcelg
  ,const word *codeg,const float4 *velrhopg,const float4 *velrhopm1g,const float3 *aceg){
  //-Reserva memoria basica.
  const unsigned n=pfin-pini;
  //-Carga posicion.
  tfloat3 *pos=new tfloat3[n];
  {
    tdouble2 *pxy=new tdouble2[n];
    double *pz=new double[n];
    cudaMemcpy(pxy,posxyg+pini,sizeof(double2)*n,cudaMemcpyDeviceToHost);
    cudaMemcpy(pz,poszg+pini,sizeof(double)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++)pos[p]=TFloat3(float(pxy[p].x),float(pxy[p].y),float(pz[p]));
    delete[] pxy;
    delete[] pz;
  }
  //-Carga idp.
  unsigned *idp=NULL;
  if(idpg){
    idp=new unsigned[n];
    cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  }
  //-Carga dcel.
  tuint3 *dcel=NULL;
  if(dcelg){
    dcel=new tuint3[n];
    unsigned *aux=new unsigned[n];
    cudaMemcpy(aux,dcelg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++)dcel[p]=TUint3(unsigned(PC__Cellx(cellcode,aux[p])),unsigned(PC__Celly(cellcode,aux[p])),unsigned(PC__Cellz(cellcode,aux[p])));
    delete[] aux;
  }
  //-Carga vel y rhop.
  tfloat3 *vel=NULL;
  float *rhop=NULL;
  if(velrhopg){
    vel=new tfloat3[n];
    rhop=new float[n];
    tfloat4 *aux=new tfloat4[n];
    cudaMemcpy(aux,velrhopg+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++){ vel[p]=TFloat3(aux[p].x,aux[p].y,aux[p].z); rhop[p]=aux[p].w; }
    delete[] aux;
  }
  //-Carga velm1 y rhopm1.
  tfloat3 *velm1=NULL;
  float *rhopm1=NULL;
  if(velrhopm1g){
    velm1=new tfloat3[n];
    rhopm1=new float[n];
    tfloat4 *aux=new tfloat4[n];
    cudaMemcpy(aux,velrhopm1g+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++){ velm1[p]=TFloat3(aux[p].x,aux[p].y,aux[p].z); rhopm1[p]=aux[p].w; }
    delete[] aux;
  }
  //-Carga ace.
  tfloat3 *ace=NULL;
  if(aceg){
    ace=new tfloat3[n];
    cudaMemcpy(ace,aceg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  }
  //-Carga type.
  byte *type=NULL;
  if(codeg){
    type=new byte[n];
    word *aux=new word[n];
    cudaMemcpy(aux,codeg+pini,sizeof(word)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++){ 
      const word cod=aux[p];
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
    delete[] aux;
  }

  //-Define campos
  JFormatFiles2::StScalarData fields[10];
  unsigned nfields=0;
  if(idp){    fields[nfields]=JFormatFiles2::DefineField("Id"    ,JFormatFiles2::UInt32  ,1,idp);    nfields++; }
  if(dcel){   fields[nfields]=JFormatFiles2::DefineField("Dcel"  ,JFormatFiles2::UInt32  ,3,dcel);   nfields++; }
  if(vel){    fields[nfields]=JFormatFiles2::DefineField("Vel"   ,JFormatFiles2::Float32 ,3,vel);    nfields++; }
  if(rhop){   fields[nfields]=JFormatFiles2::DefineField("Rhop"  ,JFormatFiles2::Float32 ,1,rhop);   nfields++; }
  if(velm1){  fields[nfields]=JFormatFiles2::DefineField("Velm1" ,JFormatFiles2::Float32 ,3,velm1);  nfields++; }
  if(rhopm1){ fields[nfields]=JFormatFiles2::DefineField("Rhopm1",JFormatFiles2::Float32 ,1,rhopm1); nfields++; }
  if(ace){    fields[nfields]=JFormatFiles2::DefineField("Ace"   ,JFormatFiles2::Float32 ,3,ace);    nfields++; }
  if(type){   fields[nfields]=JFormatFiles2::DefineField("Typex" ,JFormatFiles2::UChar8  ,1,type);   nfields++; }

  //-Genera fichero.
  JFormatFiles2::SaveVtk(fun::FileNameSec(filename,numfile),n,pos,nfields,fields);

  //-Libera memoria 
  delete[] pos;
  delete[] idp;
  delete[] dcel;
  delete[] vel;
  delete[] rhop;
  delete[] velm1;
  delete[] rhopm1;
  delete[] ace;
  delete[] type;
}

