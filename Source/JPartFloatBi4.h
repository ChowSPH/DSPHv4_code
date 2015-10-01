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
//# Clase JPartFloatBi4Save para grabar la informacion de los floatings.
//# Clase JPartFloatBi4Load para recuperar la informacion de los floatings.
//# Algunas de sus funcionalidades son:
//# - Graba datos de cabecera y de estado de floatings por PART.
//# - Comprueba cabecera y recupera datos de estado de floatings por PART.
//#
//# Cambios:
//# =========
//# - Implementacion. (04-12-2014)
//# - Implementacion independiente de StFloatingData en Types.h. (11-05-2015)
//# - Remplaza long long por llong. (01-10-2015)
//#############################################################################

/// \file JPartFloatBi4.h \brief Declares the class \ref JPartFloatBi4Save and class \ref JPartFloatBi4Load.

#ifndef _JPartFloatBi4_
#define _JPartFloatBi4_

#include "JObject.h"
#include "TypesDef.h"
#include "JBinaryData.h"
#include <string>
#include <vector>
#include <fstream>


//##############################################################################
//# JPartFloatBi4Save
//##############################################################################
class JPartFloatBi4Save : protected JObject
{
 private:
  JBinaryData *Data;      ///<Almacena la informacion general de los datos (constante para cada PART).
  JBinaryData *Part;      ///<Pertenece a Data y almacena informacion de un part (incluyendo datos de floatings).

  //-Variables de gestion.
  static const unsigned FormatVerDef=141204;    ///<Version de formato by default.
  unsigned FormatVer;    ///<Version de formato.

  bool InitialSaved;     ///<Indica si se grabo la informacion de cabecera.

  std::string AppName;   ///<Nombre de aplicacion.
  std::string Dir;       ///<Directorio de datos.
  unsigned FtCount;      ///<Numero de floatings.

  //-Datos constantes de floatings (head).
  word *HeadMkbound;
  unsigned *HeadBegin;
  unsigned *HeadCount;
  float *HeadMass;
  float *HeadRadius;

  //-Datos variables de floatings (PARTs).
  tdouble3 *PartCenter;
  tfloat3 *PartFvel;
  tfloat3 *PartFomega;
  
  unsigned Cpart;    ///<Numero de PART.

  void ResizeFtData(unsigned ftcount);
  void ClearPartData();
  static std::string GetNamePart(unsigned cpart);

 public:
  JPartFloatBi4Save();
  ~JPartFloatBi4Save();
  void Reset();
  void ResetData();
  void ResetPart();

  llong GetAllocMemory()const;
  static std::string GetFileNamePart();

  //Grabacion de datos:
  //====================
  //-Configuracion de objeto.
  void Config(std::string appname,const std::string &dir,unsigned ftcount);
  void AddHeadData(unsigned cf,word mkbound,unsigned begin,unsigned count,float mass,float radius);
  void SaveInitial();

  ////-Configuracion de parts.
  void AddPartData(unsigned cf,const tdouble3 &center,const tfloat3 &fvel,const tfloat3 &fomega);
  JBinaryData* AddPartFloat(unsigned cpart,double timestep,double demdtforce);

  ////-Grabacion de fichero.
  void SavePartFloat();
  void SavePartFloat(unsigned cpart,double timestep,double demdtforce){   AddPartFloat(cpart,timestep,demdtforce); SavePartFloat();   }
};


//##############################################################################
//# JPartFloatBi4Load
//##############################################################################
class JPartFloatBi4Load : protected JObject
{
 private:
  JBinaryData *Data;      ///<Almacena la informacion general de los datos (constante para cada PART).
  unsigned FtCount;       ///<Numero de floatings.
  unsigned PartCount;     ///<Numero de PARTs.
  JBinaryData *Part;      ///<Pertenece a Data y almacena informacion de un part (incluyendo datos de floatings).

  //-Datos constantes de floatings (head).
  word *HeadMkbound;
  unsigned *HeadBegin;
  unsigned *HeadCount;
  float *HeadMass;
  float *HeadRadius;

  //-Informacion de PART.
  double TimeStep;
  double DemDtForce;
  //-Datos variables de floatings (PARTs).
  tdouble3 *PartCenter;
  tfloat3 *PartFvel;
  tfloat3 *PartFomega;

  JBinaryDataArray* CheckArray(JBinaryData *bd,const std::string &name,JBinaryDataDef::TpData type);
  void ResetPart();
  void ResizeFtData(unsigned ftcount);
  void CheckPart()const;
  void CheckFloating(unsigned cf)const;

 public:
  JPartFloatBi4Load();
  ~JPartFloatBi4Load();
  void Reset();
  static std::string GetFileNamePart();

  void LoadFile(const std::string &dir);
  void CheckHeadData(unsigned cf,word mkbound,unsigned begin,unsigned count,float mass);

  unsigned GetFtCount()const{ return(FtCount); }
  unsigned GetCount()const{ return(PartCount); }

  word     GetHeadMkbound(unsigned cf)const{ CheckFloating(cf); return(HeadMkbound[cf]); }
  unsigned GetHeadBegin  (unsigned cf)const{ CheckFloating(cf); return(HeadBegin[cf]);   }
  unsigned GetHeadCount  (unsigned cf)const{ CheckFloating(cf); return(HeadCount[cf]);   }
  float    GetHeadMass   (unsigned cf)const{ CheckFloating(cf); return(HeadMass[cf]);    }
  float    GetHeadRadius (unsigned cf)const{ CheckFloating(cf); return(HeadRadius[cf]);  }

  void LoadPart(unsigned cpart);

  double GetPartTimeStep()const{ CheckPart(); return(TimeStep); }
  double GetPartDemDtForce()const{ CheckPart(); return(DemDtForce); }
  
  tdouble3 GetPartCenter(unsigned cf)const{ CheckFloating(cf); return(PartCenter[cf]); }
  tfloat3 GetPartFvel   (unsigned cf)const{ CheckFloating(cf); return(PartFvel[cf]);   }
  tfloat3 GetPartFomega (unsigned cf)const{ CheckFloating(cf); return(PartFomega[cf]); }
};



#endif

