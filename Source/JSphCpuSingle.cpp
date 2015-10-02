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

#include "JSphCpuSingle.h"
#include "JCellDivCpuSingle.h"
#include "JArraysCpu.h"
#include "Functions.h"
#include "JXml.h"
#include "JSphMotion.h"
#include "JPartsLoad4.h"
#include "JSphVisco.h"
#include "JWaveGen.h"

#include <climits>

using namespace std;
//==============================================================================
// Constructor.
//==============================================================================
JSphCpuSingle::JSphCpuSingle():JSphCpu(false){
  ClassName="JSphCpuSingle";
  CellDivSingle=NULL;
  PartsLoaded=NULL;
}

//==============================================================================
// Destructor.
//==============================================================================
JSphCpuSingle::~JSphCpuSingle(){
  delete CellDivSingle; CellDivSingle=NULL;
  delete PartsLoaded;   PartsLoaded=NULL;
}

//==============================================================================
// Devuelve la memoria reservada en cpu.
//==============================================================================
llong JSphCpuSingle::GetAllocMemoryCpu()const{  
  llong s=JSphCpu::GetAllocMemoryCpu();
  //Reservada en otros objetos
  if(CellDivSingle)s+=CellDivSingle->GetAllocMemory();
  if(PartsLoaded)s+=PartsLoaded->GetAllocMemory();
  return(s);
}

//==============================================================================
// Actualiza los valores maximos de memory, particles y cells.
//==============================================================================
void JSphCpuSingle::UpdateMaxValues(){
  MaxParticles=max(MaxParticles,Np);
  if(CellDivSingle)MaxCells=max(MaxCells,CellDivSingle->GetNct());
  llong m=GetAllocMemoryCpu();
  MaxMemoryCpu=max(MaxMemoryCpu,m);
}

//==============================================================================
// Carga la configuracion de ejecucion.
//==============================================================================
void JSphCpuSingle::LoadConfig(JCfgRun *cfg){
  const char met[]="LoadConfig";
  //-Carga configuracion de OpenMP
  ConfigOmp(cfg);
  //-Carga configuracion basica general
  JSph::LoadConfig(cfg);
  //-Checks compatibility of selected options.
  if(RenCorrection && UseDEM)RunException(met,"Ren correction is not implemented with Floatings-DEM.");
  Log->Print("**Special case configuration is loaded");
}

//==============================================================================
// Carga particulas del caso a procesar.
//==============================================================================
void JSphCpuSingle::LoadCaseParticles(){
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
void JSphCpuSingle::ConfigDomain(){
  const char* met="ConfigDomain";
  //-Calcula numero de particulas.
  Np=PartsLoaded->GetCount(); Npb=CaseNpb; NpbOk=Npb;
  //-Reserva memoria fija para moving y floating.
  AllocCpuMemoryFixed();
  //-Reserva memoria en Cpu para particulas.
  AllocCpuMemoryParticles(Np,0);

  //-Copia datos de particulas.
  ReserveBasicArraysCpu();
  memcpy(Posc,PartsLoaded->GetPos(),sizeof(tdouble3)*Np);
  memcpy(Idpc,PartsLoaded->GetIdp(),sizeof(unsigned)*Np);
  memcpy(Velrhopc,PartsLoaded->GetVelRhop(),sizeof(tfloat4)*Np);

  //-Calcula radio de floatings.
  if(CaseNfloat && PeriActive!=0)CalcFloatingRadius(Np,Posc,Idpc);

  //-Carga code de particulas.
  LoadCodeParticles(Np,Idpc,Codec);

  //-Libera memoria de PartsLoaded.
  delete PartsLoaded; PartsLoaded=NULL;
  //-Aplica configuracion de CellOrder.
  ConfigCellOrder(CellOrder,Np,Posc,Velrhopc);

  //-Configura division celdas.
  ConfigCellDivision();
  //-Establece dominio de simulacion local dentro de Map_Cells y calcula DomCellCode.
  SelecDomain(TUint3(0,0,0),Map_Cells);
  //-Calcula celda inicial de particulas y comprueba si hay excluidas inesperadas.
  LoadDcellParticles(Np,Codec,Posc,Dcellc);

  //-Crea objeto para divide en Gpu y selecciona un cellmode valido.
  CellDivSingle=new JCellDivCpuSingle(Stable,FtCount!=0,PeriActive,CellOrder,CellMode,Scell,Map_PosMin,Map_PosMax,Map_Cells,CaseNbound,CaseNfixed,CaseNpb,Log,DirOut);
  CellDivSingle->DefineDomain(DomCellCode,DomCelIni,DomCelFin,DomPosMin,DomPosMax);
  ConfigCellDiv((JCellDivCpu*)CellDivSingle);

  ConfigSaveData(0,1,"");

  //-Reordena particulas por celda.
  BoundChanged=true;
  RunCellDivide(true);
}

//==============================================================================
// Redimensiona el espacio reservado para particulas en CPU midiendo el
// tiempo consumido con TMC_SuResizeNp. Al terminar actualiza el divide.
//==============================================================================
void JSphCpuSingle::ResizeParticlesSize(unsigned newsize,float oversize,bool updatedivide){
  TmcStart(Timers,TMC_SuResizeNp);
  newsize+=(oversize>0? unsigned(oversize*newsize): 0);
  ResizeCpuMemoryParticles(newsize);
  TmcStop(Timers,TMC_SuResizeNp);
  if(updatedivide)RunCellDivide(true);
}

//==============================================================================
// Crea lista de nuevas particulas periodicas a duplicar.
// Con stable activado reordena lista de periodicas.
//==============================================================================
unsigned JSphCpuSingle::PeriodicMakeList(unsigned n,unsigned pini,bool stable,unsigned nmax,tdouble3 perinc,const tdouble3 *pos,const word *code,unsigned *listp)const{
  unsigned count=0;
  if(n){
    //-Inicializa tamaño de lista lspg a cero.
    listp[nmax]=0;
    for(unsigned p=0;p<n;p++){
      const unsigned p2=p+pini;
      //-Se queda con particulas normales o periodicas.
      if(CODE_GetSpecialValue(code[p2])<=CODE_PERIODIC){
        //-Obtiene posicion de particula.
        const tdouble3 ps=pos[p2];
        tdouble3 ps2=ps+perinc;
        if(Map_PosMin<=ps2 && ps2<Map_PosMax){
          unsigned cp=listp[nmax]; listp[nmax]++; if(cp<nmax)listp[cp]=p2;
        }
        ps2=ps-perinc;
        if(Map_PosMin<=ps2 && ps2<Map_PosMax){
          unsigned cp=listp[nmax]; listp[nmax]++; if(cp<nmax)listp[cp]=(p2|0x80000000);
        }
      }
    }
    count=listp[nmax];
    //-Reordena lista si es valida y stable esta activado.
    if(stable && count && count<=nmax){
      //-No hace falta porque de momento no se crea la lista usando OpenMP.
    }
  }
  return(count);
}

//==============================================================================
// Duplica la posicion de la particula indicada aplicandole un desplazamiento.
// Las particulas duplicadas se considera que siempre son validas y estan dentro
// del dominio.
// Este kernel vale para single-cpu y multi-cpu porque los calculos se hacen 
// a partir de domposmin.
// Se controla que las coordendas de celda no sobrepasen el maximo.
//==============================================================================
void JSphCpuSingle::PeriodicDuplicatePos(unsigned pnew,unsigned pcopy,bool inverse,double dx,double dy,double dz,tuint3 cellmax,tdouble3 *pos,unsigned *dcell)const{
  //-Obtiene pos de particula a duplicar.
  tdouble3 ps=pos[pcopy];
  //-Aplica desplazamiento.
  ps.x+=(inverse? -dx: dx);
  ps.y+=(inverse? -dy: dy);
  ps.z+=(inverse? -dz: dz);
  //-Calcula coordendas de celda dentro de dominio.
  unsigned cx=unsigned((ps.x-DomPosMin.x)/Scell);
  unsigned cy=unsigned((ps.y-DomPosMin.y)/Scell);
  unsigned cz=unsigned((ps.z-DomPosMin.z)/Scell);
  //-Ajusta las coordendas de celda si sobrepasan el maximo.
  cx=(cx<=cellmax.x? cx: cellmax.x);
  cy=(cy<=cellmax.y? cy: cellmax.y);
  cz=(cz<=cellmax.z? cz: cellmax.z);
  //-Graba posicion y celda de nuevas particulas.
  pos[pnew]=ps;
  dcell[pnew]=PC__Cell(DomCellCode,cx,cy,cz);
}

//==============================================================================
// Crea particulas periodicas a partir de una lista con las particulas a duplicar.
// Se presupone que todas las particulas son validas.
// Este kernel vale para single-cpu y multi-cpu porque usa domposmin. 
//==============================================================================
void JSphCpuSingle::PeriodicDuplicateVerlet(unsigned n,unsigned pini,tuint3 cellmax,tdouble3 perinc,const unsigned *listp
  ,unsigned *idp,word *code,unsigned *dcell,tdouble3 *pos,tfloat4 *velrhop,tsymatrix3f *spstau,tfloat4 *velrhopm1)const
{
  for(unsigned p=0;p<n;p++){
    const unsigned pnew=p+pini;
    const unsigned rp=listp[p];
    const unsigned pcopy=(rp&0x7FFFFFFF);
    //-Ajusta posicion y celda de nueva particula.
    PeriodicDuplicatePos(pnew,pcopy,(rp>=0x80000000),perinc.x,perinc.y,perinc.z,cellmax,pos,dcell);
    //-Copia el resto de datos.
    idp[pnew]=idp[pcopy];
    code[pnew]=CODE_SetPeriodic(code[pcopy]);
    velrhop[pnew]=velrhop[pcopy];
    velrhopm1[pnew]=velrhopm1[pcopy];
    if(spstau)spstau[pnew]=spstau[pcopy];
  }
}

//==============================================================================
// Crea particulas periodicas a partir de una lista con las particulas a duplicar.
// Se presupone que todas las particulas son validas.
// Este kernel vale para single-cpu y multi-cpu porque usa domposmin. 
//==============================================================================
void JSphCpuSingle::PeriodicDuplicateSymplectic(unsigned n,unsigned pini,tuint3 cellmax,tdouble3 perinc,const unsigned *listp
  ,unsigned *idp,word *code,unsigned *dcell,tdouble3 *pos,tfloat4 *velrhop,tsymatrix3f *spstau,tdouble3 *pospre,tfloat4 *velrhoppre)const
{
  for(unsigned p=0;p<n;p++){
    const unsigned pnew=p+pini;
    const unsigned rp=listp[p];
    const unsigned pcopy=(rp&0x7FFFFFFF);
    //-Ajusta posicion y celda de nueva particula.
    PeriodicDuplicatePos(pnew,pcopy,(rp>=0x80000000),perinc.x,perinc.y,perinc.z,cellmax,pos,dcell);
    //-Copia el resto de datos.
    idp[pnew]=idp[pcopy];
    code[pnew]=CODE_SetPeriodic(code[pcopy]);
    velrhop[pnew]=velrhop[pcopy];
    if(pospre)pospre[pnew]=pospre[pcopy];
    if(velrhoppre)velrhoppre[pnew]=velrhoppre[pcopy];
    if(spstau)spstau[pnew]=spstau[pcopy];
  }
}

//==============================================================================
// Crea particulas duplicadas de condiciones periodicas.
// Crea nuevas particulas periodicas y marca las viejas para ignorarlas.
// Las nuevas periodicas se situan a partir del Np de entrada, primero las NpbPer
// de contorno y despues las NpfPer fluidas. El Np de salida contiene tambien las
// nuevas periodicas.
//==============================================================================
void JSphCpuSingle::RunPeriodic(){
  const char met[]="RunPeriodic";
  TmcStart(Timers,TMC_SuPeriodic);
  //-Guarda numero de periodicas actuales.
  NpfPerM1=NpfPer;
  NpbPerM1=NpbPer;
  //-Marca periodicas actuales para ignorar.
  for(unsigned p=0;p<Np;p++){
    const word rcode=Codec[p];
    if(CODE_GetSpecialValue(rcode)==CODE_PERIODIC)Codec[p]=CODE_SetOutIgnore(rcode);
  }
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
          unsigned* listp=ArraysCpu->ReserveUint();
          unsigned nmax=CpuParticlesSize-1; //-Numero maximo de particulas que caben en la lista.
          //-Genera lista de nuevas periodicas.
          if(Np>=0x80000000)RunException(met,"The number of particles is too big.");//-Pq el ultimo bit se usa para marcar el sentido en que se crea la nueva periodica.
          unsigned count=PeriodicMakeList(num2,pini2,Stable,nmax,perinc,Posc,Codec,listp);
          //-Redimensiona memoria para particulas si no hay espacio suficiente y repite el proceso de busqueda.
          if(count>nmax || count+Np>CpuParticlesSize){
            ArraysCpu->Free(listp); listp=NULL;
            TmcStop(Timers,TMC_SuPeriodic);
            ResizeParticlesSize(Np+count,PERIODIC_OVERMEMORYNP,false);
            TmcStart(Timers,TMC_SuPeriodic);
          }
          else{
            run=false;
            //-Crea nuevas particulas periodicas duplicando las particulas de la lista.
            if(TStep==STEP_Verlet)PeriodicDuplicateVerlet(count,Np,DomCells,perinc,listp,Idpc,Codec,Dcellc,Posc,Velrhopc,SpsTauc,VelrhopM1c);
            if(TStep==STEP_Symplectic){
              if((PosPrec || VelrhopPrec) && (!PosPrec || !VelrhopPrec))RunException(met,"Symplectic data is invalid.") ;
              PeriodicDuplicateSymplectic(count,Np,DomCells,perinc,listp,Idpc,Codec,Dcellc,Posc,Velrhopc,SpsTauc,PosPrec,VelrhopPrec);
            }

            //-Libera lista y actualiza numero de particulas.
            ArraysCpu->Free(listp); listp=NULL;
            Np+=count;
            //-Actualiza numero de periodicas nuevas.
            if(!ctype)NpbPer+=count;
            else NpfPer+=count;
          }
        }
      }
    }
  }
  TmcStop(Timers,TMC_SuPeriodic);
}

//==============================================================================
// Ejecuta divide de particulas en celdas.
//==============================================================================
void JSphCpuSingle::RunCellDivide(bool updateperiodic){
  const char met[]="RunCellDivide";
  //-Crea nuevas particulas periodicas y marca las viejas para ignorarlas.
  if(updateperiodic && PeriActive)RunPeriodic();

  //-Inicia Divide.
  CellDivSingle->Divide(Npb,Np-Npb-NpbPer-NpfPer,NpbPer,NpfPer,BoundChanged,Dcellc,Codec,Idpc,Posc,Timers);

  //-Ordena datos de particulas
  TmcStart(Timers,TMC_NlSortData);
  CellDivSingle->SortArray(Idpc);
  CellDivSingle->SortArray(Codec);
  CellDivSingle->SortArray(Dcellc);
  CellDivSingle->SortArray(Posc);
  CellDivSingle->SortArray(Velrhopc);
  if(TStep==STEP_Verlet){
    CellDivSingle->SortArray(VelrhopM1c);
  }
  else if(TStep==STEP_Symplectic && (PosPrec || VelrhopPrec)){//En realidad solo es necesario en el divide del corrector, no en el predictor???
    if(!PosPrec || !VelrhopPrec)RunException(met,"Symplectic data is invalid.") ;
    CellDivSingle->SortArray(PosPrec);
    CellDivSingle->SortArray(VelrhopPrec);
  }
  if(TVisco==VISCO_LaminarSPS)CellDivSingle->SortArray(SpsTauc);

  //-Recupera datos del divide.
  Np=CellDivSingle->GetNpFinal();
  Npb=CellDivSingle->GetNpbFinal();
  NpbOk=Npb-CellDivSingle->GetNpbIgnore();
  //-Recupera posiciones de floatings.
  if(CaseNfloat)CalcRidp(PeriActive!=0,Np-Npb,Npb,CaseNpb,CaseNpb+CaseNfloat,Codec,Idpc,FtRidp);
  TmcStop(Timers,TMC_NlSortData);

  //-Gestion de particulas excluidas (solo fluid pq si alguna bound es excluida se genera excepcion en Divide()).
  TmcStart(Timers,TMC_NlOutCheck);
  unsigned npfout=CellDivSingle->GetNpOut();
  if(npfout){
    unsigned* idp=ArraysCpu->ReserveUint();
    tdouble3* pos=ArraysCpu->ReserveDouble3();
    tfloat3* vel=ArraysCpu->ReserveFloat3();
    float* rhop=ArraysCpu->ReserveFloat();
    unsigned num=GetParticlesData(npfout,Np,true,false,idp,pos,vel,rhop,NULL);
    AddParticlesOut(npfout,idp,pos,vel,rhop,CellDivSingle->GetNpfOutRhop(),CellDivSingle->GetNpfOutMove());
    ArraysCpu->Free(idp);
    ArraysCpu->Free(pos);
    ArraysCpu->Free(vel);
    ArraysCpu->Free(rhop);
  }
  TmcStop(Timers,TMC_NlOutCheck);
  BoundChanged=false;
}

//------------------------------------------------------------------------------
// Devuelve limites de celdas para interaccion.
//------------------------------------------------------------------------------
void JSphCpuSingle::GetInteractionCells(unsigned rcell
  ,int hdiv,const tint4 &nc,const tint3 &cellzero
  ,int &cxini,int &cxfin,int &yini,int &yfin,int &zini,int &zfin)const
{
  //-Obtiene limites de interaccion
  const int cx=PC__Cellx(DomCellCode,rcell)-cellzero.x;
  const int cy=PC__Celly(DomCellCode,rcell)-cellzero.y;
  const int cz=PC__Cellz(DomCellCode,rcell)-cellzero.z;
  //-Codigo para hdiv 1 o 2 pero no cero.
  cxini=cx-min(cx,hdiv);
  cxfin=cx+min(nc.x-cx-1,hdiv)+1;
  yini=cy-min(cy,hdiv);
  yfin=cy+min(nc.y-cy-1,hdiv)+1;
  zini=cz-min(cz,hdiv);
  zfin=cz+min(nc.z-cz-1,hdiv)+1;
}

//==============================================================================
// Aplica correccion de Ren a la presion y densidad del contorno.
//==============================================================================
void JSphCpuSingle::RunRenCorrection(){
  //-Calcula presion en contorno a partir de fluido.
  float *presskf=ArraysCpu->ReserveFloat();
  Interaction_Ren(NpbOk,CellDivSingle->GetNcells(),CellDivSingle->GetBeginCell()
    ,CellDivSingle->GetCellDomainMin(),Dcellc,Posc,PsPosc,Velrhopc,Idpc,Codec,Pressc,presskf);
  //-Recalcula valores de presion y densidad en contorno segun RenBeta.
  ComputeRenPress(NpbOk,RenCorrection,presskf,Velrhopc,Pressc);
  ArraysCpu->Free(presskf); presskf=NULL;
}

//==============================================================================
// Interaccion para el calculo de fuerzas.
//==============================================================================
void JSphCpuSingle::Interaction_Forces(TpInter tinter){
  const char met[]="Interaction_Forces";
  PreInteraction_Forces(tinter);
  TmcStart(Timers,TMC_CfForces);
  if(RenCorrection)RunRenCorrection();

  //-Interaccion Fluid-Fluid/Bound & Bound-Fluid (forces and DEM).
  float viscdt=0;
  if(Psimple)JSphCpu::InteractionSimple_Forces(Np,Npb,NpbOk,CellDivSingle->GetNcells(),CellDivSingle->GetBeginCell(),CellDivSingle->GetCellDomainMin(),Dcellc,PsPosc,Velrhopc,Idpc,Codec,Pressc,viscdt,Arc,Acec,Deltac,SpsTauc,SpsGradvelc,ShiftPosc,ShiftDetectc);
  else JSphCpu::Interaction_Forces(Np,Npb,NpbOk,CellDivSingle->GetNcells(),CellDivSingle->GetBeginCell(),CellDivSingle->GetCellDomainMin(),Dcellc,Posc,Velrhopc,Idpc,Codec,Pressc,viscdt,Arc,Acec,Deltac,SpsTauc,SpsGradvelc,ShiftPosc,ShiftDetectc);

  //-Para simulaciones 2D anula siempre la 2º componente
  if(Simulate2D)for(unsigned p=Npb;p<Np;p++)Acec[p].y=0;

  //-Añade correccion de Delta-SPH a Arg[].
  if(Deltac){
    const int ini=int(Npb),fin=int(Np),npf=int(Np-Npb);
    #ifdef _WITHOMP
      #pragma omp parallel for schedule (static) if(npf>LIMIT_COMPUTELIGHT_OMP)
    #endif
    for(int p=ini;p<fin;p++)if(Deltac[p]!=FLT_MAX)Arc[p]+=Deltac[p];
  }

  //-Calculates maximum value of ViscDt.
  ViscDtMax=viscdt;
  //-Calculates maximum value of Ace.
  AceMax=ComputeAceMax();

  TmcStop(Timers,TMC_CfForces);
}

//==============================================================================
// Devuelve valor maximo de (ace.x^2 + ace.y^2 + ace.z^2) a partir de Acec[].
// The use of OpenMP here is not efficient.
//==============================================================================
double JSphCpuSingle::ComputeAceMax(){
  float acemax=0;
  const int ini=int(Npb),fin=int(Np),npf=int(Np-Npb);
  if(!PeriActive){//-Sin condiciones periodicas.
    for(int p=ini;p<fin;p++){
      const float ace=Acec[p].x*Acec[p].x+Acec[p].y*Acec[p].y+Acec[p].z*Acec[p].z;
      acemax=max(acemax,ace);
    }
  }
  else{//-Con condiciones periodicas ignora las particulas periodicas.
    for(int p=ini;p<fin;p++)if(CODE_GetSpecialValue(Codec[p])==CODE_NORMAL){
      const float ace=Acec[p].x*Acec[p].x+Acec[p].y*Acec[p].y+Acec[p].z*Acec[p].z;
      acemax=max(acemax,ace);
    }
  }
  return(sqrt(double(acemax)));
}

//==============================================================================
// Realiza interaccion y actualizacion de particulas segun las fuerzas 
// calculadas en la interaccion usando Verlet.
//==============================================================================
double JSphCpuSingle::ComputeStep_Ver(){
  Interaction_Forces(INTER_Forces);    //-Interaccion
  const double dt=DtVariable(true);    //-Calcula nuevo dt
  DemDtForce=dt;                       //(DEM)
  if(TShifting)RunShifting(dt);        //-Shifting
  ComputeVerlet(dt);                   //-Actualiza particulas usando Verlet
  if(CaseNfloat)RunFloating(dt,false); //-Gestion de floating bodies
  PosInteraction_Forces();             //-Libera memoria de interaccion
  return(dt);
}

//==============================================================================
// Realiza interaccion y actualizacion de particulas segun las fuerzas 
// calculadas en la interaccion usando Symplectic.
//==============================================================================
double JSphCpuSingle::ComputeStep_Sym(){
  const double dt=DtPre;
  //-Predictor
  //-----------
  DemDtForce=dt*0.5f;                     //(DEM)
  Interaction_Forces(INTER_Forces);       //-Interaccion
  const double ddt_p=DtVariable(false);   //-Calcula dt del predictor
  if(TShifting)RunShifting(dt*.5);        //-Shifting
  ComputeSymplecticPre(dt);               //-Aplica Symplectic-Predictor a las particulas
  if(CaseNfloat)RunFloating(dt*.5,true);  //-Gestion de floating bodies
  PosInteraction_Forces();                //-Libera memoria de interaccion
  //-Corrector
  //-----------
  DemDtForce=dt;                          //(DEM)
  RunCellDivide(true);
  Interaction_Forces(INTER_ForcesCorr);   //Interaccion
  const double ddt_c=DtVariable(true);    //-Calcula dt del corrector
  if(TShifting)RunShifting(dt);           //-Shifting
  ComputeSymplecticCorr(dt);              //-Aplica Symplectic-Corrector a las particulas
  if(CaseNfloat)RunFloating(dt,false);    //-Gestion de floating bodies
  PosInteraction_Forces();                //-Libera memoria de interaccion

  DtPre=min(ddt_p,ddt_c);                 //-Calcula el dt para el siguiente ComputeStep
  return(dt);
}

//==============================================================================
// Calcula distancia entre pariculas floatin y centro segun condiciones periodicas.
//==============================================================================
tfloat3 JSphCpuSingle::FtPeriodicDist(const tdouble3 &pos,const tdouble3 &center,float radius)const{
  tdouble3 distd=(pos-center);
  if(PeriX && fabs(distd.x)>radius){
    if(distd.x>0)distd=distd+PeriXinc;
    else distd=distd-PeriXinc;
  }
  if(PeriY && fabs(distd.y)>radius){
    if(distd.y>0)distd=distd+PeriYinc;
    else distd=distd-PeriYinc;
  }
  if(PeriZ && fabs(distd.z)>radius){
    if(distd.z>0)distd=distd+PeriZinc;
    else distd=distd-PeriZinc;
  }
  return(ToTFloat3(distd));
}

//==============================================================================
// Calcula fuerzas sobre floatings.
//==============================================================================
void JSphCpuSingle::FtCalcForces(StFtoForces *ftoforces)const{
  const int ftcount=int(FtCount);
  #ifdef _WITHOMP
    #pragma omp parallel for schedule (guided)
  #endif
  for(int cf=0;cf<ftcount;cf++){
    const StFloatingData fobj=FtObjs[cf];
    const unsigned fpini=fobj.begin-CaseNpb;
    const unsigned fpfin=fpini+fobj.count;
    const float fradius=fobj.radius;
    const tdouble3 fcenter=fobj.center;
    const float fmassp=fobj.massp;
    //-Computes traslational and rotational velocities.
    tfloat3 face=TFloat3(0);
    tfloat3 fomegavel=TFloat3(0);
    tmatrix3f inert=TMatrix3f(0,0,0,0,0,0,0,0,0);
    //-Calcula sumatorios: face, fomegavel y inert.
    for(unsigned fp=fpini;fp<fpfin;fp++){
      int p=FtRidp[fp];
      //-Ace is initialised with the value of the gravity for all particles.
      float acex=Acec[p].x-Gravity.x,acey=Acec[p].y-Gravity.y,acez=Acec[p].z-Gravity.z;
      face.x+=acex; face.y+=acey; face.z+=acez;
      tfloat3 dist=(PeriActive? FtPeriodicDist(Posc[p],fcenter,fradius): ToTFloat3(Posc[p]-fcenter)); 
      fomegavel.x+= acez*dist.y - acey*dist.z;
      fomegavel.y+= acex*dist.z - acez*dist.x;
      fomegavel.z+= acey*dist.x - acex*dist.y;
      //inertia tensor
      inert.a11+=(float)  (dist.y*dist.y+dist.z*dist.z)*fmassp;
      inert.a12+=(float) -(dist.x*dist.y)*fmassp;
      inert.a13+=(float) -(dist.x*dist.z)*fmassp;
      inert.a21+=(float) -(dist.x*dist.y)*fmassp;
      inert.a22+=(float)  (dist.x*dist.x+dist.z*dist.z)*fmassp;
      inert.a23+=(float) -(dist.y*dist.z)*fmassp;
      inert.a31+=(float) -(dist.x*dist.z)*fmassp;
      inert.a32+=(float) -(dist.y*dist.z)*fmassp;
      inert.a33+=(float)  (dist.x*dist.x+dist.y*dist.y)*fmassp;
    }
    //-Calculates the inverse of the intertia matrix to compute the I^-1 * L= W
    tmatrix3f invinert=TMatrix3f(0,0,0,0,0,0,0,0,0);
    const float detiner=(inert.a11*inert.a22*inert.a33+inert.a12*inert.a23*inert.a31+inert.a21*inert.a32*inert.a13-(inert.a31*inert.a22*inert.a13+inert.a21*inert.a12*inert.a33+inert.a23*inert.a32*inert.a11));
    if(detiner){
      invinert.a11= (inert.a22*inert.a33-inert.a23*inert.a32)/detiner;
      invinert.a12=-(inert.a12*inert.a33-inert.a13*inert.a32)/detiner;
      invinert.a13= (inert.a12*inert.a23-inert.a13*inert.a22)/detiner;
      invinert.a21=-(inert.a21*inert.a33-inert.a23*inert.a31)/detiner;
      invinert.a22= (inert.a11*inert.a33-inert.a13*inert.a31)/detiner;
      invinert.a23=-(inert.a11*inert.a23-inert.a13*inert.a21)/detiner;
      invinert.a31= (inert.a21*inert.a32-inert.a22*inert.a31)/detiner;
      invinert.a32=-(inert.a11*inert.a32-inert.a12*inert.a31)/detiner;
      invinert.a33= (inert.a11*inert.a22-inert.a12*inert.a21)/detiner;
    }
    //-Calcula omega a partir de fomegavel y invinert.
    {
      tfloat3 omega;
      omega.x=(fomegavel.x*invinert.a11+fomegavel.y*invinert.a12+fomegavel.z*invinert.a13);
      omega.y=(fomegavel.x*invinert.a21+fomegavel.y*invinert.a22+fomegavel.z*invinert.a23);
      omega.z=(fomegavel.x*invinert.a31+fomegavel.y*invinert.a32+fomegavel.z*invinert.a33);
      fomegavel=omega;
    }
    //-Guarda resultados en ftoforces[].
    ftoforces[cf].face=face;
    ftoforces[cf].fomegavel=fomegavel;
  }
}

//==============================================================================
// Procesa floating objects.
//==============================================================================
void JSphCpuSingle::RunFloating(double dt,bool predictor){
  const char met[]="RunFloating";
  if(TimeStep>=FtPause){//-Se usa >= pq si FtPause es cero en symplectic-predictor no entraria.
    TmcStart(Timers,TMC_SuFloating);
    //-Calcula fuerzas sobre floatings.
    FtCalcForces(FtoForces);

    //-Aplica movimiento sobre floatings.
    const int ftcount=int(FtCount);
    #ifdef _WITHOMP
      #pragma omp parallel for schedule (guided)
    #endif
    for(int cf=0;cf<ftcount;cf++){
      //-Obtiene datos de floating.
      const StFloatingData fobj=FtObjs[cf];
      //-Calculo de face.
      const float fmass=fobj.mass;
      tfloat3 face=FtoForces[cf].face;
      face.x=(face.x+fmass*Gravity.x)/fmass;
      face.y=(face.y+fmass*Gravity.y)/fmass;
      face.z=(face.z+fmass*Gravity.z)/fmass;
      //-Calculo de fomega.
      tfloat3 fomega=fobj.fomega;
      {
        const tfloat3 omega=FtoForces[cf].fomegavel;
        fomega.x=float(dt*omega.x+fomega.x);
        fomega.y=float(dt*omega.y+fomega.y);
        fomega.z=float(dt*omega.z+fomega.z);
      }
      tfloat3 fvel=fobj.fvel;
      //-Anula componentes para 2D.
      if(Simulate2D){ face.y=0; fomega.x=0; fomega.z=0; fvel.y=0; }
      //-Calculo de fcenter.
      tdouble3 fcenter=fobj.center;
      fcenter.x+=dt*fvel.x;
      fcenter.y+=dt*fvel.y;
      fcenter.z+=dt*fvel.z;
      //-Calculo de fvel.
      fvel.x=float(dt*face.x+fvel.x);
      fvel.y=float(dt*face.y+fvel.y);
      fvel.z=float(dt*face.z+fvel.z);

      //-Updates floating particles.
      const float fradius=fobj.radius;
      const unsigned fpini=fobj.begin-CaseNpb;
      const unsigned fpfin=fpini+fobj.count;
      for(unsigned fp=fpini;fp<fpfin;fp++){
        const int p=FtRidp[fp];
        if(p!=UINT_MAX){
          tfloat4 *velrhop=Velrhopc+p;
          //-Calcula y graba desplazamiento de posicion.
          const double dx=dt*double(velrhop->x);
          const double dy=dt*double(velrhop->y);
          const double dz=dt*double(velrhop->z);
          UpdatePos(Posc[p],dx,dy,dz,false,p,Posc,Dcellc,Codec);
          //-Calcula y graba nueva velocidad.
          tfloat3 dist=(PeriActive? FtPeriodicDist(Posc[p],fcenter,fradius): ToTFloat3(Posc[p]-fcenter)); 
          velrhop->x=fvel.x+(fomega.y*dist.z-fomega.z*dist.y);
          velrhop->y=fvel.y+(fomega.z*dist.x-fomega.x*dist.z);
          velrhop->z=fvel.z+(fomega.x*dist.y-fomega.y*dist.x);
        }
      }

      //-Stores floating data.
      if(!predictor){
        const tdouble3 centerold=FtObjs[cf].center;
        FtObjs[cf].center=(PeriActive? UpdatePeriodicPos(fcenter): fcenter);
        FtObjs[cf].fvel=fvel;
        FtObjs[cf].fomega=fomega;
      }
    }
    TmcStop(Timers,TMC_SuFloating);
  }
}

//==============================================================================
// Inicia proceso de simulacion.
//==============================================================================
void JSphCpuSingle::Run(std::string appname,JCfgRun *cfg,JLog2 *log){
  const char* met="Run";
  if(!cfg||!log)return;
  AppName=appname; Log=log;

  //-Configura timers
  //-------------------
  TmcCreation(Timers,cfg->SvTimers);
  TmcStart(Timers,TMC_Init);
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
  ConfigRunMode(cfg);

  //-Inicializacion de variables de ejecucion
  //-------------------------------------------
  InitRun();
  UpdateMaxValues();
  PrintAllocMemory(GetAllocMemoryCpu());
  SaveData(); 
  TmcResetValues(Timers);
  TmcStop(Timers,TMC_Init);
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
    //if(Nstep>=3)break;
  }
  TimerSim.Stop(); TimerTot.Stop();

  //-Fin de simulacion
  //--------------------
  FinishRun(partoutstop);
}

//==============================================================================
// Genera los ficheros de salida de datos
//==============================================================================
void JSphCpuSingle::SaveData(){
  const bool save=(SvData!=SDAT_None&&SvData!=SDAT_Info);
  const unsigned npsave=Np-NpbPer-NpfPer; //-Resta las periodicas si las hubiera.
  TmcStart(Timers,TMC_SuSavePart);
  //-Recupera datos de particulas en orden original.
  unsigned *idp=NULL;
  tdouble3 *pos=NULL;
  tfloat3 *vel=NULL;
  float *rhop=NULL;
  if(save){
    //-Asigna memoria y recupera datos de las particulas.
    idp=ArraysCpu->ReserveUint();
    pos=ArraysCpu->ReserveDouble3();
    vel=ArraysCpu->ReserveFloat3();
    rhop=ArraysCpu->ReserveFloat();
    unsigned npnormal=GetParticlesData(Np,0,true,PeriActive!=0,idp,pos,vel,rhop,NULL);
    if(npnormal!=npsave)RunException("SaveData","The number of particles is invalid.");
  }
  //-Reune informacion adicional.
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
    infoplus.gpudata=false;
    TimerSim.Stop();
    infoplus.timesim=TimerSim.GetElapsedTimeD()/1000.;
  }
  //-Graba datos de particulas.
  const tdouble3 vdom[2]={OrderDecode(CellDivSingle->GetDomainLimits(true)),OrderDecode(CellDivSingle->GetDomainLimits(false))};
  JSph::SaveData(npsave,idp,pos,vel,rhop,1,vdom,&infoplus);
  //-Libera memoria para datos de particulas.
  ArraysCpu->Free(idp);
  ArraysCpu->Free(pos);
  ArraysCpu->Free(vel);
  ArraysCpu->Free(rhop);
  //-Graba informacion de ejecucion.
  if(TimersStep)TimersStep->SaveData();
  TmcStop(Timers,TMC_SuSavePart);
}

//==============================================================================
// Muestra y graba resumen final de ejecucion.
//==============================================================================
void JSphCpuSingle::FinishRun(bool stop){
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

