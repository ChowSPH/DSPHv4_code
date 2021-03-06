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

#ifndef _JCellDivCpu_
#define _JCellDivCpu_

//#############################################################################
//# Cambios:
//# =========
//#############################################################################

#include "Types.h"
#include "JObject.h"
#include "JSphTimersCpu.h"
#include "JLog2.h"
#include <cmath>
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>

//#define DBG_JCellDivCpu 1 //pdte

class JCellDivCpu : protected JObject
{
protected:
  const bool Stable;
  const bool Floating;
  const byte PeriActive;
  const TpCellOrder CellOrder;
  const TpCellMode CellMode;    //-Modo de division en celdas.
  const unsigned Hdiv;          //-Valor por el que se divide a DosH
  const float Scell,OvScell;
  const tdouble3 Map_PosMin,Map_PosMax,Map_PosDif;
  const tuint3 Map_Cells;
  const unsigned CaseNbound,CaseNfixed,CaseNpb;
  JLog2 *Log;
  std::string DirOut;

  bool AllocFullNct; //-Reserva memoria para el numero maximo de celdas del dominio (DomCells).
  float OverMemoryNp;//-Porcentaje que se a�ade a la reserva de memoria de Np. (def=0).
  word OverMemoryCells;//-Numero celdas que se incrementa en cada dimension reservar memoria. (def=0).

  //-Vars del dominio definido.
  unsigned DomCellCode;   //-Clave para la codificacion de la celda de posicion.
  tuint3 DomCelIni,DomCelFin;
  tdouble3 DomPosMin,DomPosMax;
  tuint3 DomCells;

  //-Memoria reservada en funcion de particulas.
  unsigned SizeNp;
  unsigned *CellPart;
  unsigned *SortPart;

  //-Memoria reservada en funcion de celdas.
  unsigned SizeNct;
  unsigned *PartsInCell;
  unsigned *BeginCell;  //-Contiene el principio de cada celda. 
  // BeginCell=[BoundOk(nct),BoundIgnore(1),Fluid(nct),BoundOut(1),FluidOut(1),BoundOutIgnore(1),FluidOutIgnore(1),END)]

  //-Variables para reordenar particulas
  byte *VSort;//-Memoria para reordenar particulas. [sizeof(tdouble3)*Np]
  word *VSortWord;//-Para ordenar vectores word (apunta a VSort).
  int *VSortInt;//-Para ordenar vectores int (apunta a VSort).
  float *VSortFloat;//-Para ordenar vectores float (apunta a VSort).
  tfloat3 *VSortFloat3;//-Para ordenar vectores tfloat3 (apunta a VSort).
  tfloat4 *VSortFloat4;//-Para ordenar vectores tfloat4 (apunta a VSort).
  tdouble3 *VSortDouble3;//-Para ordenar vectores tdouble3 (apunta a VSort).
  tsymatrix3f *VSortSymmatrix3f;//-Para ordenar vectores tsymatrix3f (apunta a VSort).

  long long MemAllocNp;  //-Mermoria reservada para particulas.
  long long MemAllocNct; //-Mermoria reservada para celdas.

  unsigned Ndiv,NdivFull;

  //-Numero de particulas por tipo al iniciar el divide.
  unsigned Npb1;
  unsigned Npf1;
  unsigned Npb2;
  unsigned Npf2;

  unsigned Nptot;  //-Numero total de particulas incluidas las que se excluyeron al terminar el divide.
  unsigned NpbOut,NpfOut,NpbOutIgnore,NpfOutIgnore;
  
  unsigned NpFinal,NpbFinal;
  unsigned NpfOutRhop,NpfOutMove,NpbIgnore;

  tuint3 CellDomainMin,CellDomainMax; //-Limites del dominio en celdas dentro de DomCells.
  unsigned Ncx,Ncy,Ncz,Nsheet,Nct;
  ullong Nctt; //-Numero total de celdas incluyendo las especiales Nctt=SizeBeginCell()
  unsigned BoxIgnore,BoxFluid,BoxBoundOut,BoxFluidOut,BoxBoundOutIgnore,BoxFluidOutIgnore;

  bool BoundLimitOk;  //-Indica que los limites del contorno ya estan calculados en BoundLimitCellMin y BoundLimitCellMax.
  tuint3 BoundLimitCellMin,BoundLimitCellMax;

  bool BoundDivideOk;   //-Indica que los limites del contorno utilizados en el divide previo fueron BoundDivideCellMin y BoundDivideCellMax.
  tuint3 BoundDivideCellMin,BoundDivideCellMax;

  bool DivideFull;  //-Indica que el divide se aplico a fluido y contorno (no solo al fluido).

  void Reset();

  //-Gestion de reserva dinamica de memoria.
  void FreeMemoryNct();
  void FreeMemoryNp();
  void FreeMemoryAll();
  void SetMemoryVSort(byte *vsort);
  void AllocMemoryNp(ullong np);
  void AllocMemoryNct(ullong nct);
  void CheckMemoryNp(unsigned npmin);
  void CheckMemoryNct(unsigned nctmin);

  ullong SizeBeginCell(ullong nct)const{ return((nct*2)+5+1); } //-[BoundOk(nct),BoundIgnore(1),Fluid(nct),BoundOut(1),FluidOut(1),BoundOutIgnore(1),FluidOutIgnore(1),END(1)]

  ullong GetAllocMemoryNp()const{ return(MemAllocNp); };
  ullong GetAllocMemoryNct()const{ return(MemAllocNct); };
  ullong GetAllocMemory()const{ return(GetAllocMemoryNp()+GetAllocMemoryNct()); };

  void VisuBoundaryOut(unsigned p,unsigned id,tdouble3 pos,word code)const;
  //tuint3 GetMapCell(const tfloat3 &pos)const;
  void LimitsCellBound(unsigned n,unsigned pini,const unsigned* dcellc,const word* codec,const unsigned* idpc,const tdouble3* posc,tuint3 &cellmin,tuint3 &cellmax)const;
  void CalcCellDomainBound(unsigned n,unsigned pini,unsigned n2,unsigned pini2,const unsigned* dcellc,const word* codec,const unsigned* idpc,const tdouble3* posc,tuint3 &cellmin,tuint3 &cellmax);
  void LimitsCellFluid(unsigned n,unsigned pini,const unsigned* dcellc,const word* codec,const unsigned* idpc,const tdouble3* posc,tuint3 &cellmin,tuint3 &cellmax,unsigned &npfoutrhop,unsigned &npfoutmove)const;
  void CalcCellDomainFluid(unsigned n,unsigned pini,unsigned n2,unsigned pini2,const unsigned* dcellc,const word* codec,const unsigned* idpc,const tdouble3* posc,tuint3 &cellmin,tuint3 &cellmax);

  unsigned CellSize(unsigned box)const{ return(BeginCell[box+1]-BeginCell[box]); }

public:
  JCellDivCpu(bool stable,bool floating,byte periactive,TpCellOrder cellorder,TpCellMode cellmode,float scell,tdouble3 mapposmin,tdouble3 mapposmax,tuint3 mapcells,unsigned casenbound,unsigned casenfixed,unsigned casenpb,JLog2 *log,std::string dirout,bool allocfullnct=true,float overmemorynp=CELLDIV_OVERMEMORYNP,word overmemorycells=CELLDIV_OVERMEMORYCELLS);
  ~JCellDivCpu();

  void DefineDomain(unsigned cellcode,tuint3 domcelini,tuint3 domcelfin,tdouble3 domposmin,tdouble3 domposmax);

  void SortArray(word *vec);
  void SortArray(unsigned *vec);
  void SortArray(float *vec);
  void SortArray(tdouble3 *vec);
  void SortArray(tfloat3 *vec);
  void SortArray(tfloat4 *vec);
  void SortArray(tsymatrix3f *vec);

  TpCellMode GetCellMode()const{ return(CellMode); }
  unsigned GetHdiv()const{ return(Hdiv); }
  float GetScell()const{ return(Scell); }

  unsigned GetNct()const{ return(Nct); }
  unsigned GetNcx()const{ return(Ncx); }
  unsigned GetNcy()const{ return(Ncy); }
  unsigned GetNcz()const{ return(Ncz); }
  tuint3 GetNcells()const{ return(TUint3(Ncx,Ncy,Ncz)); }
  unsigned GetBoxFluid()const{ return(BoxFluid); }

  tuint3 GetCellDomainMin()const{ return(CellDomainMin); }
  tuint3 GetCellDomainMax()const{ return(CellDomainMax); }
  tdouble3 GetDomainLimits(bool limitmin,unsigned slicecellmin=0)const;

  unsigned GetNpFinal()const{ return(NpFinal); }
  unsigned GetNpbFinal()const{ return(NpbFinal); }
  unsigned GetNpbIgnore()const{ return(NpbIgnore); }
  unsigned GetNpOut()const{ return(NpbOut+NpfOut); }
  unsigned GetNpbOutIgnore()const{ return(NpbOutIgnore); }
  unsigned GetNpfOutIgnore()const{ return(NpfOutIgnore); }

  unsigned GetNpfOutPos()const{ return(NpfOut-(NpfOutMove+NpfOutRhop)); }
  unsigned GetNpfOutMove()const{ return(NpfOutMove); }
  unsigned GetNpfOutRhop()const{ return(NpfOutRhop); }

  //const unsigned* GetCellPart()const{ return(CellPart); }
  const unsigned* GetBeginCell(){ return(BeginCell); }

  //bool CellNoEmpty(unsigned box,byte kind)const;
  //unsigned CellBegin(unsigned box,byte kind)const;
  //unsigned CellSize(unsigned box,byte kind)const;
};

#endif


