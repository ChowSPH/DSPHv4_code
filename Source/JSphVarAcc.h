/*
<DUALSPHYSICS>  Copyright (C) 2013 by Jose M. Dominguez, Dr Alejandro Crespo, Prof. M. Gomez Gesteira, Anxo Barreiro, Ricardo Canelas
                                      Dr Benedict Rogers, Dr Stephen Longshaw, Dr Renato Vacondio

EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

This file is part of DualSPHysics. 

DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphVarAcc.h \brief Declares the class \ref JSphVarAcc.

#ifndef _JSphVarAcc_
#define _JSphVarAcc_

#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

#include "JObject.h"
#include "Types.h"
#include "JTimer.h"
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>


//##############################################################################
//# JSphVarAccFile
//##############################################################################
/// \brief Provides the force to be applied to different blocks of particles that is loaded from files.

class JSphVarAccFile : protected JObject
{
protected:
  static const unsigned SIZEMAX=100000;
  static const unsigned SIZEINITIAL=100;

  unsigned MkValue;          ///<The MK values stored in the acceleration input file.
  tfloat3 AccCoG;            ///<The centre of gravity that will be used for angular acceleration calculations.

  unsigned AccSize;          ///<Number of acceleration values that were allocated.
  unsigned AccCount;         ///<Number of acceleration values in each input file(s).
  float *AccTime;            ///<Variable acceleration time evolution as detailed in the input file.
  tfloat3 *AccLin;           ///<Linear acceleration variable to store values as they are read from the input files.
  tfloat3 *AccAng;           ///<Angular acceleration variable to store values as they are read from the input files.

  unsigned AccIndex;         ///<Current index for variable acceleration interpolation.

  tfloat3 CurrAccLin;        ///<The current interpolated values for linear acceleration.
  tfloat3 CurrAccAng;        ///<The current interpolated values for angular acceleration.

  void Resize(unsigned size);

public:
  JSphVarAccFile();
  ~JSphVarAccFile();
  void Reset();
  long long GetAllocMemory()const;
  void LoadFile(std::string file,float tmax);
  void GetAccValues(float timestep,unsigned &mkvalue,tfloat3 &acclin,tfloat3 &accang,tfloat3 &centre);
};

//##############################################################################
//# JSphVarAcc
//##############################################################################
/// \brief Manages the application of external forces to different blocks of particles (with the same MK).

class JSphVarAcc : protected JObject
{
protected:
  std::string BaseFile;
  std::vector<JSphVarAccFile*> Files;
  long long MemSize;

public:
  JSphVarAcc();
  ~JSphVarAcc();
  void Reset();
  long long GetAllocMemory()const{ return(MemSize); }

  void Config(std::string basefile,unsigned files,float tmax);
  void GetAccValues(unsigned cfile,float timestep,unsigned &mkvalue,tfloat3 &acclin,tfloat3 &accang,tfloat3 &centre);

  std::string GetBaseFile()const{ return(BaseFile); };
  unsigned GetCount()const{ return(unsigned(Files.size())); };
};

#endif


