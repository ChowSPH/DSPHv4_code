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

#include "JCellDivCpu.h"
#include "Functions.h"
#include "JFormatFiles2.h"
#include <cfloat>
#include <climits>

using namespace std;

//==============================================================================
// Constructor.
//==============================================================================
JCellDivCpu::JCellDivCpu(bool stable,bool floating,byte periactive,TpCellOrder cellorder,TpCellMode cellmode,float scell,tdouble3 mapposmin,tdouble3 mapposmax,tuint3 mapcells,unsigned casenbound,unsigned casenfixed,unsigned casenpb,JLog2 *log,std::string dirout,bool allocfullnct,float overmemorynp,word overmemorycells):Stable(stable),Floating(floating),PeriActive(periactive),CellOrder(cellorder),CellMode(cellmode),Hdiv(cellmode==CELLMODE_2H? 1: (cellmode==CELLMODE_H? 2: 0)),Scell(scell),OvScell(1.f/scell),Map_PosMin(mapposmin),Map_PosMax(mapposmax),Map_PosDif(mapposmax-mapposmin),Map_Cells(mapcells),CaseNbound(casenbound),CaseNfixed(casenfixed),CaseNpb(casenpb),Log(log),DirOut(dirout),AllocFullNct(allocfullnct),OverMemoryNp(overmemorynp),OverMemoryCells(overmemorycells)
{
  ClassName="JCellDivCpu";
  CellPart=NULL;    SortPart=NULL;
  PartsInCell=NULL; BeginCell=NULL;
  VSort=NULL;
  Reset();
}

//==============================================================================
// Destructor.
//==============================================================================
JCellDivCpu::~JCellDivCpu(){
  Log->Printf("---> DivideFull:%u/%u",NdivFull,Ndiv);
  Reset();
}

//==============================================================================
// Initialization of variables.
//==============================================================================
void JCellDivCpu::Reset(){
  SizeNp=SizeNct=0;
  FreeMemoryAll();
  Ndiv=NdivFull=0;
  Nptot=Npb1=Npf1=Npb2=Npf2=0;
  MemAllocNp=MemAllocNct=0;
  NpbOut=NpfOut=NpbOutIgnore=NpfOutIgnore=0;
  NpFinal=NpbFinal=0;
  NpfOutRhop=NpfOutMove=NpbIgnore=0;
  CellDomainMin=TUint3(1);
  CellDomainMax=TUint3(0);
  Ncx=Ncy=Ncz=Nsheet=Nct=0;
  Nctt=0;
  BoundLimitOk=BoundDivideOk=false;
  BoundLimitCellMin=BoundLimitCellMax=TUint3(0);
  BoundDivideCellMin=BoundDivideCellMax=TUint3(0);
  DivideFull=false;
}

//==============================================================================
// Libera memoria reservada para celdas.
//==============================================================================
void JCellDivCpu::FreeMemoryNct(){
  delete[] PartsInCell;   PartsInCell=NULL;
  delete[] BeginCell;     BeginCell=NULL; 
  MemAllocNct=0;
  BoundDivideOk=false;
}

//==============================================================================
// Libera memoria reservada para celdas.
//==============================================================================
void JCellDivCpu::FreeMemoryNp(){
  delete[] CellPart;    CellPart=NULL;
  delete[] SortPart;    SortPart=NULL;
  delete[] VSort;       SetMemoryVSort(NULL);
  MemAllocNp=0;
  BoundDivideOk=false;
}

//==============================================================================
// Libera memoria reservada para particulas y celdas.
//==============================================================================
void JCellDivCpu::FreeMemoryAll(){
  FreeMemoryNct();
  FreeMemoryNp();
}

//==============================================================================
// Ajusta buffers para reordenar datos de particulas.
//==============================================================================
void JCellDivCpu::SetMemoryVSort(byte *vsort){
  VSort=vsort;
  VSortInt=(int*)VSort;        VSortWord=(word*)VSort;
  VSortFloat=(float*)VSort;    VSortFloat3=(tfloat3*)VSort;
  VSortFloat4=(tfloat4*)VSort; VSortDouble3=(tdouble3*)VSort;
  VSortSymmatrix3f=(tsymatrix3f*)VSort;
}

//==============================================================================
// Asigna memoria segun numero de particulas. 
//==============================================================================
void JCellDivCpu::AllocMemoryNp(ullong np){
  const char met[]="AllocMemoryNp";
  FreeMemoryNp();
  SizeNp=unsigned(np);
  //-Comprueba numero de particulas.
  if(np!=SizeNp)RunException(met,string("Failed memory allocation for ")+fun::UlongStr(np)+" particles.");
  //-Reserva memoria para particulas.
  MemAllocNp=0;
  try{
    CellPart=new unsigned[SizeNp];                      MemAllocNp+=sizeof(unsigned)*SizeNp;
    SortPart=new unsigned[SizeNp];                      MemAllocNp+=sizeof(unsigned)*SizeNp;
    SetMemoryVSort(new byte[sizeof(tdouble3)*SizeNp]);  MemAllocNp+=sizeof(tdouble3)*SizeNp;
  }
  catch(const std::bad_alloc){
    RunException(met,fun::PrintStr("Failed CPU memory allocation of %.1f MB for %u particles.",double(MemAllocNp)/(1024*1024),SizeNp));
  }
  //-Muestra la memoria solicitada.
  Log->Printf("**CellDiv: Requested cpu memory for %u particles: %.1f MB.",SizeNp,double(MemAllocNp)/(1024*1024));
}

//==============================================================================
// Asigna memoria segun numero de celdas. 
//==============================================================================
void JCellDivCpu::AllocMemoryNct(ullong nct){
  const char met[]="AllocMemoryNct";
  FreeMemoryNct();
  SizeNct=unsigned(nct);
  //-Comprueba numero de celdas.
  if(nct!=SizeNct)RunException(met,string("Failed GPU memory allocation for ")+fun::UlongStr(nct)+" cells.");
  //-Reserva memoria para celdas.
  MemAllocNct=0;
  const unsigned nc=(unsigned)SizeBeginCell(nct);
  try{
    PartsInCell=new unsigned[nc-1];  MemAllocNct+=sizeof(unsigned)*(nc-1);
    BeginCell=new unsigned[nc];      MemAllocNct+=sizeof(unsigned)*(nc);
  }
  catch(const std::bad_alloc){
    RunException(met,fun::PrintStr("Failed CPU memory allocation of %.1f MB for %u cells.",double(MemAllocNct)/(1024*1024),SizeNct));
  }
  //-Muestra la memoria solicitada.
  Log->Printf("**CellDiv: Requested cpu memory for %u cells (CellMode=%s): %.1f MB.",SizeNct,GetNameCellMode(CellMode),double(MemAllocNct)/(1024*1024));
}

//==============================================================================
// Comprueba la reserva de memoria para el numero indicado de particulas. 
// Si no es suficiente o no hay reserva, entonces reserva la memoria requerida.
//==============================================================================
void JCellDivCpu::CheckMemoryNp(unsigned npmin){
  if(SizeNp<npmin)AllocMemoryNp(ullong(npmin)+ullong(OverMemoryNp*npmin));
  else if(!CellPart)AllocMemoryNp(SizeNp);  
}

//==============================================================================
// Comprueba la reserva de memoria para el numero indicado de celdas. 
// Si no es suficiente o no hay reserva, entonces reserva la memoria requerida.
//==============================================================================
void JCellDivCpu::CheckMemoryNct(unsigned nctmin){
  if(SizeNct<nctmin){
    unsigned OverMemoryCells=1;
    unsigned overnct=0;
    if(OverMemoryCells>0){
      ullong nct=ullong(Ncx+1)*ullong(Ncy+1)*ullong(Ncz+1);
      ullong nctt=SizeBeginCell(nct);
      if(nctt!=unsigned(nctt))RunException("CheckMemoryNct","The number of cells is too big.");
      overnct=unsigned(nct);
    }
    AllocMemoryNct(nctmin>overnct? nctmin: overnct);
  }
  else if(!BeginCell)AllocMemoryNct(SizeNct);  
}

//==============================================================================
// Define el dominio de simulacion a usar.
//==============================================================================
void JCellDivCpu::DefineDomain(unsigned cellcode,tuint3 domcelini,tuint3 domcelfin,tdouble3 domposmin,tdouble3 domposmax){
  DomCellCode=cellcode;
  DomCelIni=domcelini;
  DomCelFin=domcelfin;
  DomPosMin=domposmin;
  DomPosMax=domposmax;
  DomCells=DomCelFin-DomCelIni;
}

//==============================================================================
// Visualiza la informacion de una particula de contorno excluida.
//==============================================================================
void JCellDivCpu::VisuBoundaryOut(unsigned p,unsigned id,tdouble3 pos,word code)const{
  string info="particle boundary out> type:";
  word tp=CODE_GetType(code);
  if(tp==CODE_TYPE_FIXED)info=info+"Fixed";
  else if(tp==CODE_TYPE_MOVING)info=info+"Moving";
  else if(tp==CODE_TYPE_FLOATING)info=info+"Floating";
  info=info+" cause:";
  word out=CODE_GetSpecialValue(code);
  if(out==CODE_OUTMOVE)info=info+"Speed";
  else if(out==CODE_OUTPOS)info=info+"Position";
  else info=info+"???";
  Log->PrintDbg(info+fun::PrintStr(" p:%u id:%u pos:(%f,%f,%f)",p,id,pos.x,pos.y,pos.z));
}

//==============================================================================
// Devuelve coordenadas de celda a partir de una posicion.
//==============================================================================
//tuint3 JCellDivCpu::GetMapCell(const tfloat3 &pos)const{
//  float dx=pos.x-MapPosMin.x,dy=pos.y-MapPosMin.y,dz=pos.z-MapPosMin.z;
//  unsigned cx=unsigned(dx*OvScell),cy=unsigned(dy*OvScell),cz=unsigned(dz*OvScell);
//  return(TUint3(cx,cy,cz));
//}

//==============================================================================
// Calcula celda minima y maxima de las particulas validas.
// En code[] ya estan marcadas las particulas excluidas.
// En caso de no haber ninguna particula valida el minimo sera mayor que el maximo.
// Si encuentra alguna particula excluida genera excepcion mostrando su info.
//==============================================================================
void JCellDivCpu::LimitsCellBound(unsigned n,unsigned pini,const unsigned* dcellc,const word* codec,const unsigned* idpc,const tdouble3* posc,tuint3 &cellmin,tuint3 &cellmax)const{
  tuint3 cmin=TUint3(1);
  tuint3 cmax=TUint3(0);
  unsigned nerr=0;
  const unsigned pfin=pini+n;
  for(unsigned p=pini;p<pfin;p++){
    unsigned rcell=dcellc[p];
    const unsigned cx=PC__Cellx(DomCellCode,rcell);
    const unsigned cy=PC__Celly(DomCellCode,rcell);
    const unsigned cz=PC__Cellz(DomCellCode,rcell);
    const word rcode=codec[p];
    const word rcodsp=CODE_GetSpecialValue(rcode);
    if(rcodsp<CODE_OUTIGNORE){ //-Particula no excluida.
      if(cmin.x>cx)cmin.x=cx;
      if(cmin.y>cy)cmin.y=cy;
      if(cmin.z>cz)cmin.z=cz;
      if(cmax.x<cx)cmax.x=cx;
      if(cmax.y<cy)cmax.y=cy;
      if(cmax.z<cz)cmax.z=cz;
    }
    else if(rcodsp>CODE_OUTIGNORE){
      if(nerr<10)VisuBoundaryOut(p,idpc[p],OrderDecodeValue(CellOrder,posc[p]),rcode);
      nerr++;
    }
  }
  if(nerr)RunException("LimitsCellBound","Some boundary particle was found outside the domain.");
  cellmin=cmin;
  cellmax=cmax;
}

//==============================================================================
// Calcula posiciones minimas y maximas del rango de particulas Bound indicado.
// En code[] ya estan marcadas las particulas excluidas.
//==============================================================================
void JCellDivCpu::CalcCellDomainBound(unsigned n,unsigned pini,unsigned n2,unsigned pini2,const unsigned* dcellc,const word* codec,const unsigned* idpc,const tdouble3* posc,tuint3 &cellmin,tuint3 &cellmax){
  tuint3 cmin,cmax;
  LimitsCellBound(n,pini,dcellc,codec,idpc,posc,cmin,cmax);
  cellmin=(cmin.x>cmax.x? DomCells: cmin);
  cellmax=(cmin.x>cmax.x? TUint3(0): cmax);
  if(n2){
    LimitsCellBound(n2,pini2,dcellc,codec,idpc,posc,cmin,cmax);
    cmin=(cmin.x>cmax.x? DomCells: cmin);
    cmax=(cmin.x>cmax.x? TUint3(0): cmax);
    cellmin=MinValues(cellmin,cmin);
    cellmax=MaxValues(cellmax,cmax);
  }
}

//==============================================================================
// Calcula celda minima y maxima de las particulas validas.
// En code[] ya estan marcadas las particulas excluidas.
// En caso de no haber ninguna particula valida el minimo sera mayor que el maximo.
// Si encuentra alguna particula floating excluida genera excepcion mostrando su info.
//==============================================================================
void JCellDivCpu::LimitsCellFluid(unsigned n,unsigned pini,const unsigned* dcellc,const word* codec,const unsigned* idpc,const tdouble3* posc,tuint3 &cellmin,tuint3 &cellmax,unsigned &npfoutrhop,unsigned &npfoutmove)const{
  unsigned noutrhop=0,noutmove=0;
  tuint3 cmin=TUint3(1);
  tuint3 cmax=TUint3(0);
  unsigned nerr=0;
  const unsigned pfin=pini+n;
  for(unsigned p=pini;p<pfin;p++){
    unsigned rcell=dcellc[p];
    const unsigned cx=PC__Cellx(DomCellCode,rcell);
    const unsigned cy=PC__Celly(DomCellCode,rcell);
    const unsigned cz=PC__Cellz(DomCellCode,rcell);
    const word rcode=codec[p];
    const word rcodsp=CODE_GetSpecialValue(rcode);
    if(rcodsp<CODE_OUTIGNORE){ //-Particula no excluida.
      if(cmin.x>cx)cmin.x=cx;
      if(cmin.y>cy)cmin.y=cy;
      if(cmin.z>cz)cmin.z=cz;
      if(cmax.x<cx)cmax.x=cx;
      if(cmax.y<cy)cmax.y=cy;
      if(cmax.z<cz)cmax.z=cz;
    }
    else if(rcodsp>CODE_OUTIGNORE){
      if(Floating && CODE_GetType(rcode)==CODE_TYPE_FLOATING){
        if(nerr<10)VisuBoundaryOut(p,idpc[p],OrderDecodeValue(CellOrder,posc[p]),codec[p]);
        nerr++;
      }
      if(rcodsp==CODE_OUTRHOP)noutrhop++;
      else if(rcodsp==CODE_OUTMOVE)noutmove++;
    }
  }
  if(nerr)RunException("LimitsCellFluid","Some floating particle was found outside the domain.");
  cellmin=cmin;
  cellmax=cmax;
  npfoutrhop+=noutrhop;
  npfoutmove+=noutmove;
}

//==============================================================================
// Calcula posiciones minimas y maximas del rango de particulas Fluid indicado.
// Ignora particulas excluidas que ya estan marcadas en code[].
//==============================================================================
void JCellDivCpu::CalcCellDomainFluid(unsigned n,unsigned pini,unsigned n2,unsigned pini2,const unsigned* dcellc,const word* codec,const unsigned* idpc,const tdouble3* posc,tuint3 &cellmin,tuint3 &cellmax){
  tuint3 cmin,cmax;
  LimitsCellFluid(n,pini,dcellc,codec,idpc,posc,cmin,cmax,NpfOutRhop,NpfOutMove);
  cellmin=(cmin.x>cmax.x? DomCells: cmin);
  cellmax=(cmin.x>cmax.x? TUint3(0): cmax);
  if(n2){
    LimitsCellFluid(n2,pini2,dcellc,codec,idpc,posc,cmin,cmax,NpfOutRhop,NpfOutMove);
    cmin=(cmin.x>cmax.x? DomCells: cmin);
    cmax=(cmin.x>cmax.x? TUint3(0): cmax);
    cellmin=MinValues(cellmin,cmin);
    cellmax=MaxValues(cellmax,cmax);
  }
  //char cad[512]; sprintf(cad,"CalcDomainFluid> cell:(%s)-(%s)",fun::Uint3Str(cellmin).c_str(),fun::Uint3Str(cellmax).c_str()); Log->PrintDbg(cad);
}

//==============================================================================
// Reordena datos de todas las particulas.
//==============================================================================
void JCellDivCpu::SortArray(word *vec){
  if(DivideFull){
    for(unsigned p=0;p<Nptot;p++)VSortWord[p]=vec[SortPart[p]];
    memcpy(vec,VSortWord,sizeof(word)*Nptot);
  }
  else{
    for(unsigned p=NpbFinal;p<Nptot;p++)VSortWord[p]=vec[SortPart[p]];
    memcpy(vec+NpbFinal,VSortWord+NpbFinal,sizeof(word)*(Nptot-NpbFinal));
  }
}
//==============================================================================
void JCellDivCpu::SortArray(unsigned *vec){
  if(DivideFull){
    for(unsigned p=0;p<Nptot;p++)VSortInt[p]=vec[SortPart[p]];
    memcpy(vec,VSortInt,sizeof(unsigned)*Nptot);
  }
  else{
    for(unsigned p=NpbFinal;p<Nptot;p++)VSortInt[p]=vec[SortPart[p]];
    memcpy(vec+NpbFinal,VSortInt+NpbFinal,sizeof(unsigned)*(Nptot-NpbFinal));
  }
}
//==============================================================================
void JCellDivCpu::SortArray(float *vec){
  if(DivideFull){
    for(unsigned p=0;p<Nptot;p++)VSortFloat[p]=vec[SortPart[p]];
    memcpy(vec,VSortFloat,sizeof(float)*Nptot);
  }
  else{
    for(unsigned p=NpbFinal;p<Nptot;p++)VSortFloat[p]=vec[SortPart[p]];
    memcpy(vec+NpbFinal,VSortFloat+NpbFinal,sizeof(float)*(Nptot-NpbFinal));
  }
}
//==============================================================================
void JCellDivCpu::SortArray(tdouble3 *vec){
  if(DivideFull){
    for(unsigned p=0;p<Nptot;p++)VSortDouble3[p]=vec[SortPart[p]];
    memcpy(vec,VSortDouble3,sizeof(tdouble3)*Nptot);
  }
  else{
    for(unsigned p=NpbFinal;p<Nptot;p++)VSortDouble3[p]=vec[SortPart[p]];
    memcpy(vec+NpbFinal,VSortDouble3+NpbFinal,sizeof(tdouble3)*(Nptot-NpbFinal));
  }
}
//==============================================================================
void JCellDivCpu::SortArray(tfloat3 *vec){
  if(DivideFull){
    for(unsigned p=0;p<Nptot;p++)VSortFloat3[p]=vec[SortPart[p]];
    memcpy(vec,VSortFloat3,sizeof(tfloat3)*Nptot);
  }
  else{
    for(unsigned p=NpbFinal;p<Nptot;p++)VSortFloat3[p]=vec[SortPart[p]];
    memcpy(vec+NpbFinal,VSortFloat3+NpbFinal,sizeof(tfloat3)*(Nptot-NpbFinal));
  }
}
//==============================================================================
void JCellDivCpu::SortArray(tfloat4 *vec){
  if(DivideFull){
    for(unsigned p=0;p<Nptot;p++)VSortFloat4[p]=vec[SortPart[p]];
    memcpy(vec,VSortFloat4,sizeof(tfloat4)*Nptot);
  }
  else{
    for(unsigned p=NpbFinal;p<Nptot;p++)VSortFloat4[p]=vec[SortPart[p]];
    memcpy(vec+NpbFinal,VSortFloat4+NpbFinal,sizeof(tfloat4)*(Nptot-NpbFinal));
  }
}
//==============================================================================
void JCellDivCpu::SortArray(tsymatrix3f *vec){
  if(DivideFull){
    for(unsigned p=0;p<Nptot;p++)VSortSymmatrix3f[p]=vec[SortPart[p]];
    memcpy(vec,VSortSymmatrix3f,sizeof(tsymatrix3f)*Nptot);
  }
  else{
    for(unsigned p=NpbFinal;p<Nptot;p++)VSortSymmatrix3f[p]=vec[SortPart[p]];
    memcpy(vec+NpbFinal,VSortSymmatrix3f+NpbFinal,sizeof(tsymatrix3f)*(Nptot-NpbFinal));
  }
}

//==============================================================================
// Devuelve limites actuales del dominio.
//==============================================================================
tdouble3 JCellDivCpu::GetDomainLimits(bool limitmin,unsigned slicecellmin)const{
  tuint3 celmin=GetCellDomainMin(),celmax=GetCellDomainMax();
  if(celmin.x>celmax.x)celmin.x=celmax.x=0; else celmax.x++;
  if(celmin.y>celmax.y)celmin.y=celmax.y=0; else celmax.y++;
  if(celmin.z>celmax.z)celmin.z=celmax.z=slicecellmin; else celmax.z++;
  double scell=double(Scell);
  tdouble3 pmin=DomPosMin+TDouble3(scell*celmin.x,scell*celmin.y,scell*celmin.z);
  tdouble3 pmax=DomPosMin+TDouble3(scell*celmax.x,scell*celmax.y,scell*celmax.z);
  return(limitmin? pmin: pmax);
}

////==============================================================================
//// Indica si la celda esta vacia o no.
////==============================================================================
//bool JCellDivCpu::CellNoEmpty(unsigned box,byte kind)const{
//#ifdef DBG_JCellDivCpu
//  if(box>=Nct)RunException("CellNoEmpty","Celda no valida.");
//#endif
//  if(kind==2)box+=BoxFluid;
//  return(BeginCell[box]<BeginCell[box+1]);
//}

////==============================================================================
//// Devuelve la primera particula de la celda solicitada.
////==============================================================================
//unsigned JCellDivCpu::CellBegin(unsigned box,byte kind)const{
//#ifdef DBG_JCellDivCpu
//  if(box>Nct)RunException("CellBegin","Celda no valida.");
//#endif
//  return(BeginCell[(kind==1? box: box+BoxFluid)]);
//}

////==============================================================================
//// Devuelve el numero de particulas de la celda solicitada.
////==============================================================================
//unsigned JCellDivCpu::CellSize(unsigned box,byte kind)const{
//#ifdef DBG_JCellDivCpu
//  if(box>Nct)RunException("CellSize","Celda no valida.");
//#endif
//  if(kind==2)box+=BoxFluid;
//  return(BeginCell[box+1]-BeginCell[box]);
//}



