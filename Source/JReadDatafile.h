/*
 <DUALSPHYSICS>  Copyright (c) 2015, Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//# Cambios:
//# =========
//# Clase para facilitar la lectura de ficheros de datos en texto. (03-11-2015)
//# - Los valores pueden estar delimitados por ',' ';' o tabs-spaces.
//# - El delimitador se identifica automaticamente.
//# - Cuenta numero de lineas de datos y de comentarios.
//# - Las lineas de comentarios empiezan por '#'.
//# - Permite extraer valores int, unsigned, float, double y string.
//# - Mantiene columna y fila del ultimo valor.
//# - En caso de usar el separador \t los espacios y tabulaciones seguidas se
//#   remplazan por una sola tabulacion. Tambien elimina espacions y tabulaciones
//#   al principio y final de la linea. (16-12-2015)
//# - Ignora lineas vacias al final del fichero. (16-12-2015)
//# - Error corregido en ProcessSpaces(). (17-12-2015)
//#############################################################################

#ifndef _JReadDatafile_
#define _JReadDatafile_

#include <string>
#include <vector>
#include "JObject.h"
#include "TypesDef.h"

//##############################################################################
//# JReadDatafile
//##############################################################################
class JReadDatafile  : protected JObject
{
private:
  std::string File;      ///< Name of file.

  unsigned SizeFile;     ///< Size of file.
  unsigned Size;         ///< Size of data.
  char *Data;            ///< Data from file.

  int LineCount;         ///< Number of lines.
  unsigned *LineBegin;   ///< Inicio de cada linea [LineCount+1].
  int RemLineCount;      ///< Number of remark lines.

  std::string Sep;       ///< Value separator.

  int ReadLin;
  int ReadLinValue;
  std::string ReadLine;
  std::string ReadValue;

  void ProcessSpaces();
  void ProcessLines();
  void ResetReadLine();

public:
  JReadDatafile();
  ~JReadDatafile();
  void Reset();

  void LoadFile(const std::string &file,unsigned maxsize=1048576000);

  unsigned Lines()const{ return(LineCount); }
  unsigned RemLines()const{ return(RemLineCount); }

  void SetReadLine(int line);

  std::string GetLine(int line)const;
  std::string ReadNextValue(bool in_line=false);

  double ReadNextDouble(bool in_line=false);
  tdouble3 ReadNextDouble3(bool in_line=false);
  
  float ReadNextFloat(bool in_line=false){    return(float(ReadNextDouble(in_line)));      }
  tfloat3 ReadNextFloat3(bool in_line=false){ return(ToTFloat3(ReadNextDouble3(in_line))); }

  int ReadNextInt(bool in_line=false);
  unsigned ReadNextUnsigned(bool in_line=false){ return(unsigned(ReadNextInt(in_line))); }

  int GetReadLin()const{           return(ReadLin);      }
  int GetReadLinValue()const{      return(ReadLinValue); }
  std::string GetReadValue()const{ return(ReadValue);    }
};

#endif


