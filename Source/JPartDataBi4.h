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
//# Clase para la gestion de ficheros en formato BI4.
//# Algunas de sus funcionalidades son:
//# - Permite la lectura y escritura de particulas en formato BI4.
//# - Solo almacena las particulas validas. Las excluidas se gestionan mediante
//#   otra clase.
//# - Permite gestionar los datos en un unico fichero o en multiples piezas.
//# - Puede almacenar la posicion como tfloat3 o tdouble3.
//# - Puede almacenar el id como unsigned de 32 o 64 bits.
//# - Esta clase esta basada en el uso de JBinaryData por lo que se puede
//#   almacenar cualquier valor o array extra.
//#
//# Cambios:
//# =========
//# - Implementacion. (12/11/2013 <-> 13/11/2013)
//# - Nuevo metodo GetFileData() para detectar si los datos usan una o varias 
//#   piezas y de paso obtener el nombre del fichero. (02/12/2013)
//# - Nuevas vars NpDynamic, ReuseIds, NpTotal y IdMax para permitir numero de
//#   particulas dinamico. (26/12/2013)
//# - Ahora el metodo SaveFileInfo() graba los datos generales al principio
//#   del fichero. (12/01/2014)
//# - Grabacion de datos propios de Splitting. (26/01/2014)
//# - Remplaza long long por llong. (01-10-2015)
//#############################################################################

/// \file JPartDataBi4.h \brief Declares the class \ref JPartDataBi4.

#ifndef _JPartDataBi4_
#define _JPartDataBi4_

#include "JObject.h"
#include "TypesDef.h"
#include "JBinaryData.h"
#include <string>
#include <vector>
#include <fstream>

//##############################################################################
//# JPartDataBi4
//##############################################################################
class JPartDataBi4 : protected JObject
{
 public:
  typedef enum{ DIV_None=0,DIV_X=1,DIV_Y=2,DIV_Z=3,DIV_Unknown=99 }TpAxisDiv; 
  typedef enum{ PERI_None=0,PERI_X=1,PERI_Y=2,PERI_Z=4,PERI_XY=3,PERI_XZ=5,PERI_YZ=6,PERI_Unknown=96 }TpPeri; 

 private:
  JBinaryData *Data;      ///<Almacena la informacion general de los datos (constante para cada PART).
  JBinaryData *Part;      ///<Pertenece a Data y almacena informacion de un part (incluyendo datos de particulas).

  //-Variables de gestion.
  static const unsigned FormatVerDef=130825;    ///<Version de formato by default.
  unsigned FormatVer;        ///<Version de formato.

  std::string Dir;   ///<Directorio de datos.
  unsigned Piece;    ///<Numero de parte.
  unsigned Npiece;   ///<Numero total de partes.
  unsigned Cpart;    ///<Numero de PART.

  static std::string GetNamePart(unsigned cpart);
  void AddPartData(unsigned npok,const unsigned *idp,const ullong *idpd,const tfloat3 *pos,const tdouble3 *posd,const tfloat3 *vel,const float *rhop);
  void SaveFileData(std::string fname);
  unsigned GetPiecesFile(std::string file)const;
  void LoadFileData(std::string file,unsigned cpart,unsigned piece,unsigned npiece);

 public:
  JPartDataBi4();
  ~JPartDataBi4();
  void Reset();
  void ResetData();
  void ResetPart();

  llong GetAllocMemory()const;
  static std::string GetFileNamePart(unsigned cpart,unsigned piece=0,unsigned npiece=1);
  static std::string GetFileNameCase(const std::string &casename,unsigned piece=0,unsigned npiece=1);
  static std::string GetFileNameInfo(unsigned piece=0,unsigned npiece=1);
  static std::string GetFileData(std::string casename,std::string dirname,unsigned cpart,byte &npiece);

  //Grabacion de datos:
  //====================
  //-Configuracion de objeto.
  void ConfigBasic(unsigned piece,unsigned npiece,std::string runcode,std::string appname,bool data2d,const std::string &dir);
  void ConfigParticles(ullong casenp,ullong casenfixed,ullong casenmoving,ullong casenfloat,ullong casenfluid,tdouble3 caseposmin,tdouble3 caseposmax,bool npdynamic=false,bool reuseids=false);
  void ConfigCtes(double dp,double h,double b,double rhop0,double gamma,double massbound,double massfluid);
  void ConfigSimMap(tdouble3 mapposmin,tdouble3 mapposmax);
  void ConfigSimPeri(TpPeri periactive,tdouble3 perixinc,tdouble3 periyinc,tdouble3 perizinc);
  void ConfigSimDiv(TpAxisDiv axisdiv);
  void ConfigSplitting(bool splitting);

  //-Configuracion de parts.
  JBinaryData* AddPartInfo(unsigned cpart,double timestep,unsigned npok,unsigned nout,unsigned step,double runtime,tdouble3 domainmin,tdouble3 domainmax,ullong nptotal=0,ullong idmax=0);
  void AddPartData(unsigned npok,const unsigned *idp,const tfloat3 *pos,const tfloat3 *vel,const float *rhop){   AddPartData(npok,idp,NULL,pos,NULL,vel,rhop);   }
  void AddPartData(unsigned npok,const unsigned *idp,const tdouble3 *posd,const tfloat3 *vel,const float *rhop){ AddPartData(npok,idp,NULL,NULL,posd,vel,rhop);  }
  void AddPartData(unsigned npok,const ullong *idpd,const tfloat3 *pos,const tfloat3 *vel,const float *rhop){    AddPartData(npok,NULL,idpd,pos,NULL,vel,rhop);  }
  void AddPartData(unsigned npok,const ullong *idpd,const tdouble3 *posd,const tfloat3 *vel,const float *rhop){  AddPartData(npok,NULL,idpd,NULL,posd,vel,rhop); }
  void AddPartDataSplitting(unsigned npok,const float *splitmass,const float *splithvar);

  //-Grabacion de fichero.
  void SaveFileCase(std::string casename);
  void SaveFilePart();
  void SaveFileInfo();

  //Carga de datos:
  //================
  //-Carga de fichero.
  unsigned GetPiecesFileCase(std::string dir,std::string casename)const;
  unsigned GetPiecesFilePart(std::string dir,unsigned cpart)const;
  void LoadFileCase(std::string dir,std::string casename,unsigned piece=0,unsigned npiece=1);
  void LoadFilePart(std::string dir,unsigned cpart,unsigned piece=0,unsigned npiece=1);

  //Obtencion de datos basicos:
  //============================
  JBinaryData* GetData()const;
  unsigned GetPiece()const{ return(Piece); } 
  unsigned GetNpiece()const{ return(Npiece); } 
  std::string Get_RunCode()const{ return(GetData()->GetvText("RunCode")); } 
  std::string Get_Date()const{    return(GetData()->GetvText("Date"));    } 
  std::string Get_AppName()const{ return(GetData()->GetvText("AppName")); } 
  bool Get_Data2d()const{         return(GetData()->GetvBool("Data2d"));  } 
  bool Get_Splitting()const{      return(GetData()->GetvBool("Splitting",true,false));  } 

  ullong Get_CaseNp()const{       return(GetData()->GetvUllong("CaseNp"));      } 
  ullong Get_CaseNfixed()const{   return(GetData()->GetvUllong("CaseNfixed"));  } 
  ullong Get_CaseNmoving()const{  return(GetData()->GetvUllong("CaseNmoving")); } 
  ullong Get_CaseNfloat()const{   return(GetData()->GetvUllong("CaseNfloat"));  } 
  ullong Get_CaseNfluid()const{   return(GetData()->GetvUllong("CaseNfluid"));  } 
  tdouble3 Get_CasePosMin()const{ return(GetData()->GetvDouble3("CasePosMin")); }
  tdouble3 Get_CasePosMax()const{ return(GetData()->GetvDouble3("CasePosMax")); }
  bool Get_NpDynamic()const{      return(GetData()->GetvBool("NpDynamic",true,false));  } 
  bool Get_ReuseIds()const{       return(GetData()->GetvBool("ReuseIds",true,false));   } 

  double Get_Dp()const{           return(GetData()->GetvDouble("Dp"));         }
  double Get_H()const{            return(GetData()->GetvDouble("H"));          }
  double Get_B()const{            return(GetData()->GetvDouble("B"));          }
  double Get_Rhop0()const{        return(GetData()->GetvDouble("Rhop0"));      }
  double Get_Gamma()const{        return(GetData()->GetvDouble("Gamma"));      }
  double Get_MassBound()const{    return(GetData()->GetvDouble("MassBound"));  }
  double Get_MassFluid()const{    return(GetData()->GetvDouble("MassFluid"));  }

  tdouble3 Get_MapPosMin()const{  return(GetData()->GetvDouble3("MapPosMin")); }
  tdouble3 Get_MapPosMax()const{  return(GetData()->GetvDouble3("MapPosMax")); }

  TpPeri Get_PeriActive()const{   return((TpPeri)GetData()->GetvInt("PeriActive")); }
  tdouble3 Get_PeriXinc()const{   return(GetData()->GetvDouble3("PeriXinc"));       }
  tdouble3 Get_PeriYinc()const{   return(GetData()->GetvDouble3("PeriYinc"));       }
  tdouble3 Get_PeriZinc()const{   return(GetData()->GetvDouble3("PeriZinc"));       }

  TpAxisDiv Get_AxisDiv()const{   return((TpAxisDiv)GetData()->GetvInt("AxisDiv")); }

  //Obtencion de datos del PART:
  //=============================
  JBinaryData* GetPart()const;
  double Get_TimeStep()const{     return(GetPart()->GetvDouble("TimeStep"));   }
  unsigned Get_Npok()const{       return(GetPart()->GetvUint("Npok"));         }
  unsigned Get_Nout()const{       return(GetPart()->GetvUint("Nout"));         }
  unsigned Get_Step()const{       return(GetPart()->GetvUint("Step"));         }
  double Get_RunTime()const{      return(GetPart()->GetvDouble("RunTime"));    }
  tdouble3 Get_DomainMin()const{  return(GetPart()->GetvDouble3("DomainMin")); }
  tdouble3 Get_DomainMax()const{  return(GetPart()->GetvDouble3("DomainMax")); }
  ullong Get_NpTotal()const{      return(GetPart()->GetvUllong("NpTotal"));    }
  ullong Get_IdMax()const{        return(GetPart()->GetvUllong("IdMax"));      }

  //Obtencion de arrays del PART:
  //==============================
  JBinaryDataArray* GetArray(std::string name)const;
  bool ArrayExists(std::string name)const;
  unsigned Get_ArrayCount(std::string name)const{ return(GetArray(name)->GetCount()); }
  bool Get_IdpSimple()const{ return(ArrayExists("Idp")); }
  bool Get_PosSimple()const{ return(ArrayExists("Pos")); }
  unsigned Get_Idp(unsigned size,unsigned *data)const{    return(GetArray("Idp")->GetDataCopy(size,data));  }
  unsigned Get_Idpd(unsigned size,ullong *data)const{     return(GetArray("Idpd")->GetDataCopy(size,data)); }
  unsigned Get_Pos(unsigned size,tfloat3 *data)const{     return(GetArray("Pos")->GetDataCopy(size,data));  }
  unsigned Get_Posd(unsigned size,tdouble3 *data)const{   return(GetArray("Posd")->GetDataCopy(size,data)); }
  unsigned Get_Vel(unsigned size,tfloat3 *data)const{     return(GetArray("Vel")->GetDataCopy(size,data));  }
  unsigned Get_Rhop(unsigned size,float *data)const{      return(GetArray("Rhop")->GetDataCopy(size,data)); }
  unsigned Get_SplitMass(unsigned size,float *data)const{ return(GetArray("SplitMass")->GetDataCopy(size,data)); }
  unsigned Get_SplitHvar(unsigned size,float *data)const{ return(GetArray("SplitHvar")->GetDataCopy(size,data)); }
};


#endif




