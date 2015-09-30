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

#include "JCellDivGpuSingle.h"
#include "JCellDivGpuSingle_ker.h"
#include "Functions.h"
#ifdef DG_JCellDivGpu
  #include "JFormatFiles2.h"
  #include "JBuffer.h"
#endif

using namespace std;

//==============================================================================
// Constructor.
//==============================================================================
JCellDivGpuSingle::JCellDivGpuSingle(bool stable,bool floating,byte periactive,TpCellOrder cellorder,TpCellMode cellmode,float scell,tdouble3 mapposmin,tdouble3 mapposmax,tuint3 mapcells,unsigned casenbound,unsigned casenfixed,unsigned casenpb,JLog2 *log,std::string dirout):JCellDivGpu(stable,floating,periactive,cellorder,cellmode,scell,mapposmin,mapposmax,mapcells,casenbound,casenfixed,casenpb,log,dirout){
  ClassName="JCellDivGpuSingle";
}

//==============================================================================
// Calcula limites del dominio en celdas ajustando al fluido (CellDomainMin/Max). 
//==============================================================================
void JCellDivGpuSingle::CalcCellDomain(const unsigned *dcellg,const word* codeg){
  //-Calcula dominio del contorno.
  tuint3 celbmin,celbmax;
  if(!BoundLimitOk){
    CalcCellDomainBound(Npb1,0,Npb2,Npb1+Npf1,dcellg,codeg,celbmin,celbmax);
    BoundLimitOk=true; BoundLimitCellMin=celbmin; BoundLimitCellMax=celbmax;
  } 
  else{ celbmin=BoundLimitCellMin; celbmax=BoundLimitCellMax; }
  //-Calcula dominio del fluido.
  tuint3 celfmin,celfmax;
  CalcCellDomainFluid(Npf1,Npb1,Npf2,Npb1+Npf1+Npb2,dcellg,codeg,celfmin,celfmax);
  //-Calcula dominio ajustando al contorno y al fluido (con halo de 2h). 
  MergeMapCellBoundFluid(celbmin,celbmax,celfmin,celfmax,CellDomainMin,CellDomainMax);
}

//==============================================================================
// Combina limite de celdas de contorno y fluido con limites de mapa.
// Con UseFluidDomain=TRUE se queda con el dominio del fluido mas 2h si hay 
// contorno, en caso contrario se queda con el dominio que incluya fluido y
// contorno.
// En caso de que el dominio sea nulo CellDomainMin=CellDomainMax=(0,0,0).
//==============================================================================
void JCellDivGpuSingle::MergeMapCellBoundFluid(const tuint3 &celbmin,const tuint3 &celbmax,const tuint3 &celfmin,const tuint3 &celfmax,tuint3 &celmin,tuint3 &celmax)const{
  //char cad[256]; sprintf(cad,"celb=(%u,%u,%u)-(%u,%u,%u)  Npb:%u",celbmin.x,celbmin.y,celbmin.z,celbmax.x,celbmax.y,celbmax.z,Npb); Log->Print(cad);
  //if(UseFluidDomain){
    celmin=TUint3(max(min(celbmin.x,celfmin.x),(celfmin.x>=Hdiv? celfmin.x-Hdiv: 0)),max(min(celbmin.y,celfmin.y),(celfmin.y>=Hdiv? celfmin.y-Hdiv: 0)),max(min(celbmin.z,celfmin.z),(celfmin.z>=Hdiv? celfmin.z-Hdiv: 0)));
    celmax=TUint3(min(max(celbmax.x,celfmax.x),celfmax.x+Hdiv),min(max(celbmax.y,celfmax.y),celfmax.y+Hdiv),min(max(celbmax.z,celfmax.z),celfmax.z+Hdiv));
 // }
 // else{
    //celmin=MinValues(celbmin,celfmin);
    //celmax=MaxValues(celbmax,celfmax);
 // }
  if(celmax.x>=DomCells.x)celmax.x=DomCells.x-1;
  if(celmax.y>=DomCells.y)celmax.y=DomCells.y-1;
  if(celmax.z>=DomCells.z)celmax.z=DomCells.z-1;
  if(celmin.x>celmax.x||celmin.y>celmax.y||celmin.z>celmax.z){ celmin=celmax=TUint3(0,0,0); }
}

//==============================================================================
// Calcula numero de celdas a partir de (CellDomainMin/Max). 
// Obtiene localizacion de celdas especiales.
//==============================================================================
void JCellDivGpuSingle::PrepareNct(){
  //-Calcula numero de celdas.
  Ncx=CellDomainMax.x-CellDomainMin.x+1;
  Ncy=CellDomainMax.y-CellDomainMin.y+1;
  Ncz=CellDomainMax.z-CellDomainMin.z+1;
  //printf("======  ncx:%u ncy:%u ncz:%u\n",Ncx,Ncy,Ncz);
  Nsheet=Ncx*Ncy; Nct=Nsheet*Ncz; Nctt=SizeBeginEndCell(Nct);
  if(Nctt!=unsigned(Nctt))RunException("PrepareNct","The number of cells is too big.");
  BoxIgnore=Nct; 
  BoxFluid=BoxIgnore+1; 
  BoxBoundOut=BoxFluid+Nct; 
  BoxFluidOut=BoxBoundOut+1; 
  BoxBoundOutIgnore=BoxFluidOut+1;
  BoxFluidOutIgnore=BoxBoundOutIgnore+1;
  //Log->Printf("--->PrepareNct> BoxIgnore:%u BoxFluid:%u BoxBoundOut:%u BoxFluidOut:%u",BoxIgnore,BoxFluid,BoxBoundOut,BoxFluidOut);
  //Log->Printf("--->PrepareNct> BoxBoundOutIgnore:%u BoxFluidOutIgnore:%u",BoxBoundOutIgnore,BoxFluidOutIgnore);
}

//==============================================================================
// Calcula celda de cada particula (CellPart[]) a partir de dcell[], todas las
// particulas excluidas ya fueron marcadas en code[].
// Asigna valores consecutivos a SortPart[].
//==============================================================================
void JCellDivGpuSingle::PreSort(const unsigned *dcellg,const word *codeg){
  if(DivideFull)cudiv::PreSortFull(Nptot,DomCellCode,dcellg,codeg,CellDomainMin,TUint3(Ncx,Ncy,Ncz),CellPart,SortPart,Log);
  else cudiv::PreSortFluid(Npf1,Npb1,DomCellCode,dcellg,codeg,CellDomainMin,TUint3(Ncx,Ncy,Ncz),CellPart,SortPart,Log);
#ifdef DG_JCellDivGpu
  tfloat4 *poscellh=new tfloat4[Np];
  unsigned *num=new unsigned[Np];
  unsigned *cellparth=new unsigned[Np];
  cudaMemcpy(poscellh,poscellg,sizeof(float4)*Np,cudaMemcpyDeviceToHost);
  for(unsigned p=0;p<Np;p++)num[p]=p;
  cudaMemcpy(cellparth,CellPart,sizeof(unsigned)*(Np),cudaMemcpyDeviceToHost);
//  unsigned nctot2=Nctot*2;
  char cad[512];
  for(unsigned p=0;p<Np;p++){
    if(cellparth[p]>=Nctt){
      sprintf(cad,"PreSort> Valor no valido de CellPart. cellpart[%u]=%u",p,cellparth[p]); Log->Print(cad);
    }  
    else if(p>=Npb&&cellparth[p]<BoxFluid){
      sprintf(cad,"PreSort> Valor de fluida no valido para CellPart. cellpart[%u]=%u",p,cellparth[p]); Log->Print(cad);
    }  
  }
  //JBuffer buf(1024*1024,1024*512);
  //buf.InStr("POINTSDATA"); buf.InUint(Np); buf.InFloat3Vec(Np,posh);
  //buf.InStr("CellPart:unsigned_int");      buf.InUintVec(Np,cellparth);
  //buf.InStr("Num:unsigned_int");           buf.InUintVec(Np,num);
  //buf.InStr("END"); 
  //JFormatFiles2::PointsToVtk(DirOut+"_CellPart.vtk",&buf);

  delete[] poscellh;
  delete[] num;
  delete[] cellparth;
#endif
}

//==============================================================================
// Inicia proceso de Divide: Calcula limites de dominio y calcula nueva posicion
// para cada particula (SortPart).
// El valor np incluye las periodicas bound y fluid (npbper y npfper).
// Las floating se tratan como si fuesen fluido (tanto al ser excluidas como ignoradas).
//==============================================================================
void JCellDivGpuSingle::Divide(unsigned npb1,unsigned npf1,unsigned npb2,unsigned npf2,bool boundchanged,const unsigned *dcellg,const word* codeg,TimersGpu timers,const double2 *posxy,const double *posz,const unsigned *idp){
  const char met[]="Divide";
  DivideFull=false;
  TmgStart(timers,TMG_NlLimits);

  //-Establece numero de particulas.
  Npb1=npb1; Npf1=npf1; Npb2=npb2; Npf2=npf2;
  Nptot=Npb1+Npf1+Npb2+Npf2;
  NpbOut=NpfOut=NpbOutIgnore=NpfOutIgnore=0;
  NpFinal=NpbFinal=0;
  NpfOutRhop=NpfOutMove=NpbIgnore=0;
  //printf("---> Npb1:%u  Npf1:%u  Npb2:%u  Npf2:%u\n",Npb1,Npf1,Npb2,Npf2);

  //-Comprueba si hay memoria reservada y si es suficiente para Nptot.
  CheckMemoryNp(Nptot);

  //-Si la posicion del contorno cambia o hay condiciones periodicas es necesario recalcular limites y reordenar todas las particulas. 
  if(boundchanged || PeriActive){
    BoundLimitOk=BoundDivideOk=false;
    BoundLimitCellMin=BoundLimitCellMax=TUint3(0);
    BoundDivideCellMin=BoundDivideCellMax=TUint3(0);
  }

  //-Calcula limites del dominio.
  CalcCellDomain(dcellg,codeg);
  //-Calcula numero de celdas para el divide y comprueba reserva de memoria para celdas.
  PrepareNct();
  //-Comprueba si hay memoria reservada y si es suficiente para Nptot.
  CheckMemoryNct(Nct);
  TmgStop(timers,TMG_NlLimits);

  //-Determina si el divide afecta a todas las particulas.
  //-BoundDivideOk se vuelve false al reservar o liberar memoria para particulas o celdas.
  if(!BoundDivideOk || BoundDivideCellMin!=CellDomainMin || BoundDivideCellMax!=CellDomainMax){
    DivideFull=true;
    BoundDivideOk=true; BoundDivideCellMin=CellDomainMin; BoundDivideCellMax=CellDomainMax;
  }
  else DivideFull=false;
//  if(DivideFull)Log->PrintDbg("--> DivideFull=TRUE"); else Log->PrintDbg("--> DivideFull=FALSE");
//  Log->PrintDbg(string("--> CellDomain:%s")+fun::Uint3RangeStr(CellDomainMin,CellDomainMax));

  //-Calcula CellPart[] y asigna valores consecutivos a SortPart[].
  TmgStart(timers,TMG_NlPreSort);
  PreSort(dcellg,codeg);
  TmgStop(timers,TMG_NlPreSort);

  //-Ordena CellPart y SortPart en funcion de la celda.
  TmgStart(timers,TMG_NlRadixSort);
  if(DivideFull)cudiv::Sort(CellPart,SortPart,Nptot,Stable);
  else cudiv::Sort(CellPart+Npb1,SortPart+Npb1,Nptot-Npb1,Stable);
  TmgStop(timers,TMG_NlRadixSort);

  //-Calcula particula inicial y final de cada celda (BeginEndCell).
  TmgStart(timers,TMG_NlCellBegin);
  cudiv::CalcBeginEndCell(DivideFull,Nptot,Npb1,unsigned(SizeBeginEndCell(Nct)),BoxFluid,CellPart,BeginEndCell);

  //-Calcula numeros de particulas.
  NpbIgnore=CellSize(BoxIgnore);
  unsigned beginendcell[8];
  CellBeginEnd(BoxBoundOut,8,beginendcell);
  NpbOut=beginendcell[1]-beginendcell[0];
  NpfOut=beginendcell[3]-beginendcell[2];
  NpbOutIgnore=beginendcell[5]-beginendcell[4];
  NpfOutIgnore=beginendcell[7]-beginendcell[6];
  //printf("---> Nct:%u  BoxBoundOut:%u  SizeBeginEndCell:%u\n",Nct,BoxBoundOut,SizeBeginEndCell(Nct));
  //printf("---> NpbIgnore:%u  NpbOut:%u  NpfOut:%u  NpfOutIgnore:%u\n",NpbIgnore,NpbOut,NpfOut,NpfOutIgnore);
  NpFinal=Nptot-NpbOut-NpfOut-NpbOutIgnore-NpfOutIgnore;
  NpbFinal=Npb1+Npb2-NpbOut-NpbOutIgnore;

  Ndiv++;
  if(DivideFull)NdivFull++;
  TmgStop(timers,TMG_NlCellBegin);
  CheckCudaError(met,"Error in NL construction.");
}



