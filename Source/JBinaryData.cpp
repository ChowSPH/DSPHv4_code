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

/// \file JBinaryData.cpp \brief Implements the class \ref JBinaryData

#include "JBinaryData.h"
#include "Functions.h"

#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <sstream>

using namespace std;

const std::string JBinaryData::CodeItemDef="\nITEM\n";
const std::string JBinaryData::CodeValuesDef="\nVALUES";
const std::string JBinaryData::CodeArrayDef="\nARRAY";

//##############################################################################
//# JBinaryDataDef
//##############################################################################
//==============================================================================
/// Devuelve tipo de datos en texto.
//==============================================================================
std::string JBinaryDataDef::TypeToStr(TpData type){
  string tx="";
  switch(type){
    case JBinaryDataDef::DatText:     tx="text";     break;
    case JBinaryDataDef::DatBool:     tx="bool";     break;
    case JBinaryDataDef::DatChar:     tx="char";     break;
    case JBinaryDataDef::DatUchar:    tx="uchar";    break;
    case JBinaryDataDef::DatShort:    tx="short";    break;
    case JBinaryDataDef::DatUshort:   tx="ushort";   break;
    case JBinaryDataDef::DatInt:      tx="int";      break;
    case JBinaryDataDef::DatUint:     tx="uint";     break;
    case JBinaryDataDef::DatLlong:    tx="llong";    break;
    case JBinaryDataDef::DatUllong:   tx="ullong";   break;
    case JBinaryDataDef::DatFloat:    tx="float";    break;
    case JBinaryDataDef::DatDouble:   tx="double";   break;
    case JBinaryDataDef::DatInt3:     tx="int3";     break;
    case JBinaryDataDef::DatUint3:    tx="uint3";    break;
    case JBinaryDataDef::DatFloat3:   tx="float3";   break;
    case JBinaryDataDef::DatDouble3:  tx="double3";  break;
  }
  return(tx);
}

//==============================================================================
/// Devuelve tama�o del tipo de datos.
//==============================================================================
size_t JBinaryDataDef::SizeOfType(TpData type){
  size_t ret=0;
  switch(type){
    //case JBinaryDataDef::DatText:
    case JBinaryDataDef::DatBool:     ret=sizeof(int);             break;
    case JBinaryDataDef::DatChar:     ret=sizeof(char);            break;
    case JBinaryDataDef::DatUchar:    ret=sizeof(unsigned char);   break;
    case JBinaryDataDef::DatShort:    ret=sizeof(short);           break;
    case JBinaryDataDef::DatUshort:   ret=sizeof(unsigned short);  break;
    case JBinaryDataDef::DatInt:      ret=sizeof(int);             break;
    case JBinaryDataDef::DatUint:     ret=sizeof(unsigned);        break;
    case JBinaryDataDef::DatLlong:    ret=sizeof(llong);           break;
    case JBinaryDataDef::DatUllong:   ret=sizeof(ullong);          break;
    case JBinaryDataDef::DatFloat:    ret=sizeof(float);           break;
    case JBinaryDataDef::DatDouble:   ret=sizeof(double);          break;
    case JBinaryDataDef::DatInt3:     ret=sizeof(tint3);           break;
    case JBinaryDataDef::DatUint3:    ret=sizeof(tuint3);          break;
    case JBinaryDataDef::DatFloat3:   ret=sizeof(tfloat3);         break;
    case JBinaryDataDef::DatDouble3:  ret=sizeof(tdouble3);        break;
  }
  return(ret);
}

//##############################################################################
//# JBinaryDataArray
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JBinaryDataArray::JBinaryDataArray(JBinaryData* parent,const std::string &name,JBinaryDataDef::TpData type):Type(type){ 
  ClassName="JBinaryDataArray";
  Parent=parent;
  Name=name;
  Hide=false;
  Pointer=NULL;
  ExternalPointer=false;
  Count=Size=0;
  ClearFileData();
}

//==============================================================================
/// Destructor.
//==============================================================================
JBinaryDataArray::~JBinaryDataArray(){
  FreeMemory();
}

//==============================================================================
/// Devuelve la cantidad de memoria reservada.
//==============================================================================
llong JBinaryDataArray::GetAllocMemory()const{
  return(Pointer&&!ExternalPointer? Size*JBinaryDataDef::SizeOfType(Type): 0);
}

//==============================================================================
/// Cambia nombre de array comprobando que no existe otro array o value con el mismo nombre.
//==============================================================================
void JBinaryDataArray::SetName(const std::string &name){
  const char met[]="SetName";
  if(Parent->ExistsValue(name))RunException(met,"There is already a value with the name given.");
  if(Parent->GetArray(name)!=NULL)RunException(met,"There is already an array with the name given.");
  if(Parent->GetItem(name)!=NULL)RunException(met,"There is already an item with the name given.");
  Name=name;
}

//==============================================================================
/// Libera memoria asignada al puntero indicado.
//==============================================================================
void JBinaryDataArray::FreePointer(void* ptr)const{
  if(ptr)switch(Type){
    case JBinaryDataDef::DatText:     delete[] (string*)ptr;          break;
    case JBinaryDataDef::DatBool:     delete[] (bool*)ptr;            break;
    case JBinaryDataDef::DatInt:      delete[] (int*)ptr;             break;
    case JBinaryDataDef::DatUint:     delete[] (unsigned int*)ptr;    break;
    case JBinaryDataDef::DatChar:     delete[] (char*)ptr;            break;
    case JBinaryDataDef::DatUchar:    delete[] (unsigned char*)ptr;   break;
    case JBinaryDataDef::DatShort:    delete[] (short*)ptr;           break;
    case JBinaryDataDef::DatUshort:   delete[] (unsigned short*)ptr;  break;
    case JBinaryDataDef::DatLlong:    delete[] (llong*)ptr;           break;
    case JBinaryDataDef::DatUllong:   delete[] (ullong*)ptr;          break;
    case JBinaryDataDef::DatFloat:    delete[] (float*)ptr;           break;
    case JBinaryDataDef::DatDouble:   delete[] (double*)ptr;          break;
    case JBinaryDataDef::DatInt3:     delete[] (tint3*)ptr;           break;  
    case JBinaryDataDef::DatUint3:    delete[] (tuint3*)ptr;          break;
    case JBinaryDataDef::DatFloat3:   delete[] (tfloat3*)ptr;         break;
    case JBinaryDataDef::DatDouble3:  delete[] (tdouble3*)ptr;        break;
    default: RunException("FreePointer","Type of array invalid.");
  }
}

//==============================================================================
/// Devuelve puntero con la memoria asiganda.
//==============================================================================
void* JBinaryDataArray::AllocPointer(unsigned size)const{
  const char met[]="AllocPointer";
  void* ptr=NULL;
  if(size){
    try{
      switch(Type){
        case JBinaryDataDef::DatText:     ptr=new string[size];          break;
        case JBinaryDataDef::DatBool:     ptr=new bool[size];            break;
        case JBinaryDataDef::DatInt:      ptr=new int[size];             break;
        case JBinaryDataDef::DatUint:     ptr=new unsigned int[size];    break;
        case JBinaryDataDef::DatChar:     ptr=new char[size];            break;
        case JBinaryDataDef::DatUchar:    ptr=new unsigned char[size];   break;
        case JBinaryDataDef::DatShort:    ptr=new short[size];           break;
        case JBinaryDataDef::DatUshort:   ptr=new unsigned short[size];  break;
        case JBinaryDataDef::DatLlong:    ptr=new llong[size];           break;
        case JBinaryDataDef::DatUllong:   ptr=new ullong[size];          break;
        case JBinaryDataDef::DatFloat:    ptr=new float[size];           break;
        case JBinaryDataDef::DatDouble:   ptr=new double[size];          break;
        case JBinaryDataDef::DatInt3:     ptr=new tint3[size];           break;
        case JBinaryDataDef::DatUint3:    ptr=new tuint3[size];          break;
        case JBinaryDataDef::DatFloat3:   ptr=new tfloat3[size];         break;
        case JBinaryDataDef::DatDouble3:  ptr=new tdouble3[size];        break;
        default: RunException(met,"Type of array invalid.");
      }
    }
    catch(const std::bad_alloc){
      RunException(met,"Cannot allocate the requested memory.");
    }
  }
  return(ptr);
}

//==============================================================================
/// Comprueba memoria disponible y redimensiona array si hace falta.
/// Si es ExternalPointer no permite redimensionar la memoria asignada.
//==============================================================================
void JBinaryDataArray::CheckMemory(unsigned count,bool resize){
  const char met[]="CheckMemory";
  if(count){
    //-Reserva memoria si fuese necesario.
    if(!Pointer){
      if(!resize)RunException(met,"Memory no allocated.");
      AllocMemory(count);
    }
    if(Count+count>Size){
      if(ExternalPointer)RunException(met,"Allocated memory in external pointer is not enough.");
      if(!resize)RunException(met,"Allocated memory is not enough.");
      AllocMemory(Count+count,true);
    }
  }
}

//==============================================================================
/// Extrae datos del ptr indicado.
//==============================================================================
void JBinaryDataArray::OutData(unsigned &count,unsigned size,const byte *ptr,byte *dat,unsigned sdat)const{
  const unsigned count2=count+sdat;
  if(count2>size)RunException("OutData","Overflow in reading data.");
  memcpy(dat,ptr+count,sdat);
  count=count2;
}

//==============================================================================
/// Extrae string de ptr.
//==============================================================================
std::string JBinaryDataArray::OutStr(unsigned &count,unsigned size,const byte *ptr)const{
  unsigned len=OutUint(count,size,ptr);
  string tex;
  tex.resize(len);
  const unsigned count2=count+len;
  if(count2>size)RunException("OutStr","Overflow in reading data.");
  memcpy((char*)tex.c_str(),ptr+count,len);
  count=count2;
  return(tex);
}


//==============================================================================
/// Libera memoria asignada.
//==============================================================================
void JBinaryDataArray::FreeMemory(){
  if(Pointer&&!ExternalPointer)FreePointer(Pointer);
  Pointer=NULL;
  ExternalPointer=false;
  Count=0;
  Size=0;
}

//==============================================================================
/// Asigna memoria para los elementos indicados.
//==============================================================================
void JBinaryDataArray::AllocMemory(unsigned size,bool savedata){
  if(Count&&savedata&&size){
    if(ExternalPointer)RunException("AllocMemory","External pointer can not be resized.");
    const unsigned count2=min(Count,size);
    void *ptr=AllocPointer(size);
    if(Type==JBinaryDataDef::DatText){//-Array de string.
      string *strings1=(string*)Pointer;
      string *strings2=(string*)ptr;
      for(unsigned c=0;c<count2;c++)strings2[c]=strings1[c];
    }
    else memcpy((byte*)ptr,(byte*)Pointer,JBinaryDataDef::SizeOfType(Type)*count2);
    FreeMemory();
    Pointer=ptr;
    Count=count2;
    Size=size;
  }
  else{
    FreeMemory();
    Size=size;
    if(Size)Pointer=AllocPointer(Size);
  }
}

//==============================================================================
/// Asigna memoria para los elementos indicados.
//==============================================================================
void JBinaryDataArray::ConfigExternalMemory(unsigned size,void* pointer){
  FreeMemory();
  ExternalPointer=true;
  Pointer=pointer;
  Size=size;
  Count=0;
}

//==============================================================================
/// Configura acceso a datos en fichero.
//==============================================================================
void JBinaryDataArray::ConfigFileData(llong filepos,unsigned datacount,unsigned datasize){
  FreeMemory();
  FileDataPos=filepos; FileDataCount=datacount; FileDataSize=datasize;
}

//==============================================================================
/// Borra datos de acceso a datos en fichero.
//==============================================================================
void JBinaryDataArray::ClearFileData(){
  FileDataPos=-1; FileDataCount=FileDataSize=0;
}

//==============================================================================
/// Carga contenido de fichero abierto con OpenFileStructure().
//==============================================================================
void JBinaryDataArray::ReadFileData(bool resize){
  const char met[]="ReadFileData";
  //printf("ReadFileData Parent_name:[%s] p:%p\n",Parent->GetName().c_str(),Parent);
  //printf("ReadFileData Parent2_name:[%s] p:%p\n",(Parent->GetParent()? Parent->GetParent()->GetName().c_str(): "none"),Parent->GetParent());
  //printf("ReadFileData root_name:[%s] p:%p\n",Parent->GetItemRoot()->GetName().c_str(),Parent->GetItemRoot());
  ifstream *pf=Parent->GetItemRoot()->GetFileStructure();
  if(!pf||!pf->is_open())RunException(met,"The file with data is not available.");
  //printf("ReadFileData[%s]> fpos:%llu count:%u size:%u\n",Name.c_str(),FileDataPos,FileDataCount,FileDataSize);
  if(FileDataPos<0)RunException(met,"The access information to data file is not available.");
  pf->seekg(FileDataPos,ios::beg);
  ReadData(FileDataCount,FileDataSize,pf,resize);
}

//==============================================================================
/// A�ade elementos al array de un fichero.
/// Si es ExternalPointer no permite redimensionar la memoria asignada.
//==============================================================================
void JBinaryDataArray::ReadData(unsigned count,unsigned size,std::ifstream *pf,bool resize){
  if(count){
    //-Reserva memoria si fuese necesario.
    CheckMemory(count,resize);
    //-Carga datos de fichero.
    if(GetType()==JBinaryDataDef::DatText){//-Array de strings.
      byte *buf=new byte[size];
      pf->read((char*)buf,size);
      unsigned cbuf=0;
      for(unsigned c=0;c<count;c++)AddText(OutStr(cbuf,size,buf),false);
      delete[] buf;
    }
    else{
      const unsigned stype=(unsigned)JBinaryDataDef::SizeOfType(Type);
      const unsigned sdat=stype*count;
      const unsigned cdat=stype*Count;
      pf->read(((char*)Pointer)+cdat,sdat);
      Count+=count;
    }
  }
}

//==============================================================================
/// A�ade elementos al array.
/// Si es ExternalPointer no permite redimensionar la memoria asignada.
//==============================================================================
void JBinaryDataArray::AddData(unsigned count,const void* data,bool resize){
  if(count){
    //-Reserva memoria si fuese necesario.
    CheckMemory(count,resize);
    //-A�ade datos al puntero.
    if(Type==JBinaryDataDef::DatText){
      string *strings=(string*)Pointer;
      string *strings2=(string*)data;
      for(unsigned c=0;c<count;c++)strings[Count+c]=strings2[c];
    }
    else{
      const unsigned stype=(unsigned)JBinaryDataDef::SizeOfType(Type);
      const unsigned sdat=stype*count;
      const unsigned cdat=stype*Count;
      memcpy(((byte*)Pointer)+cdat,(byte*)data,sdat);
    }
    Count+=count;
  }
}

//==============================================================================
/// Guarda datos como contenido del array.
//==============================================================================
void JBinaryDataArray::SetData(unsigned count,const void* data,bool externalpointer){
  FreeMemory();
  if(externalpointer){
    Pointer=(void*)data;
    ExternalPointer=true;
    Size=count;
    Count=count;
  }
  else AddData(count,data,true);
}

//==============================================================================
/// A�ade un string al array.
/// Si es ExternalPointer no permite redimensionar la memoria asignada.
//==============================================================================
void JBinaryDataArray::AddText(const std::string &str,bool resize){
  if(Type!=JBinaryDataDef::DatText)RunException("AddText","Type of array is not Text.");
  unsigned count=1;
  if(count){
    //-Reserva memoria si fuese necesario.
    CheckMemory(count,resize);
    //-A�ade string al array.
    string *strings=(string*)Pointer;
    strings[Count]=str;
    Count+=count;
  }
}

//==============================================================================
/// A�ade array de strings al array.
/// Si es ExternalPointer no permite redimensionar la memoria asignada.
//==============================================================================
void JBinaryDataArray::AddTexts(unsigned count,const std::string *strs,bool resize){
  if(Type!=JBinaryDataDef::DatText)RunException("AddTexts","Type of array is not Text.");
  if(count)AddData(count,(const void*)strs,resize);
}

//==============================================================================
/// Devuelve puntero de datos comprobando que tenga datos.
//==============================================================================
const void* JBinaryDataArray::GetDataPointer()const{
  if(!DataInPointer())RunException("GetDataPointer","There are not available data in pointer.");
  return(Pointer);
}

//==============================================================================
/// Copia datos de Pointer o FileData al puntero indicado y devuelve el numero
/// de elementos.
//==============================================================================
unsigned JBinaryDataArray::GetDataCopy(unsigned size,void* pointer)const{
  const char met[]="GetDataCopy";
  if(!DataInPointer()&&!DataInFile())RunException(met,"There are not available data in Pointer or FileData.");
  const size_t stype=JBinaryDataDef::SizeOfType(GetType());
  if(!stype)RunException(met,"Type of array is invalid for this function.");
  unsigned count=0;
  if(DataInPointer()){
    count=GetCount();
    if(size>=count)memcpy(pointer,Pointer,stype*count);
  }
  else{
    count=FileDataCount;
    if(size>=count){
      ifstream *pf=Parent->GetItemRoot()->GetFileStructure();
      if(!pf||!pf->is_open())RunException(met,"The file with data is not available.");
      pf->seekg(FileDataPos,ios::beg);
      count=FileDataCount;
      pf->read((char*)pointer,stype*count);
    }
  }
  if(size<count)RunException(met,"Size of array is not enough to store all data.");
  return(count);
}


//##############################################################################
//# JBinaryData
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JBinaryData::JBinaryData(std::string name):Name(name){ 
  ClassName="JBinaryData";
  Parent=NULL;
  FileStructure=NULL;
  ValuesData=NULL;
  ValuesCacheReset();
  HideAll=HideValues=false;
  FmtFloat="%.7E";
  FmtDouble="%.15E";
}

//==============================================================================
/// Constructor de copias.
//==============================================================================
JBinaryData::JBinaryData(const JBinaryData &src){
  ClassName="JBinaryData";
  Parent=NULL;
  FileStructure=NULL;
  ValuesData=NULL;
  ValuesCacheReset();
  *this=src;
}

//==============================================================================
/// Destructor.
//==============================================================================
JBinaryData::~JBinaryData(){
  Clear();
}

//==============================================================================
/// Sobrecarga del operador de asignacion.
//==============================================================================
JBinaryData& JBinaryData::operator=(const JBinaryData &src){
  if(this!=&src){
    unsigned size=src.GetSizeDataConst(true);
    byte *dat=new byte[size];
    src.SaveDataConst(size,dat,true);
    LoadData(size,dat);
    delete[] dat;
  }
  return(*this);
}

//==============================================================================
/// Elimina todo el contenido (values, arrays e items).
//==============================================================================
void JBinaryData::Clear(){
  RemoveValues();
  RemoveArrays();
  RemoveItems();
  ValuesCacheReset();
  CloseFileStructure();
}

//==============================================================================
/// Devuelve la cantidad de memoria reservada por el objeto y subobjetos.
//==============================================================================
llong JBinaryData::GetAllocMemory()const{
  llong s=0; 
  for(unsigned c=0;c<Arrays.size();c++)s+=Arrays[c]->GetAllocMemory(); //-Memoria de arrays no externos.
  s+=ValuesSize;                                                       //-Memoria para cache de values.
  for(unsigned c=0;c<Items.size();c++)s+=Items[c]->GetAllocMemory();   //-Memoria de Items descencientes.
  return(s);
}

//==============================================================================
/// Elimina cache de values.
//==============================================================================
void JBinaryData::ValuesCacheReset(){
  delete[] ValuesData; ValuesData=NULL;
  ValuesSize=0;
  ValuesModif=true;
}

//==============================================================================
/// Prepara cache de values.
//==============================================================================
void JBinaryData::ValuesCachePrepare(bool down){
  if(ValuesModif){
    ValuesCacheReset();
    ValuesSize=GetSizeValues();
    ValuesData=new byte[ValuesSize];
    unsigned count=0;
    SaveValues(count,ValuesSize,ValuesData);
    ValuesModif=false;
  }
  if(down)for(unsigned c=0;c<Items.size();c++)Items[c]->ValuesCachePrepare(true);
}

//==============================================================================
/// Comprueba existencia y tipo de valor. Devuelve posicion de valor (-1 no exite).
//==============================================================================
int JBinaryData::CheckGetValue(const std::string &name,bool optional,JBinaryDataDef::TpData type)const{
  const char met[]="CheckGetValue";
  int idx=GetValueIndex(name);
  //if(idx<0&&GetArrayIndex(name)>=0)RunException(met,string("The value ")+name+" is an array.");
  if(!optional&&idx<0)RunException(met,string("Value ")+name+" not found.");
  if(idx>=0&&Values[idx].type!=type)RunException(met,string("Type of value ")+name+" invalid.");
  return(idx);
}

//==============================================================================
/// Comprueba tipo de valor, y sino existe lo crea. Devuelve posicion de valor.
//==============================================================================
int JBinaryData::CheckSetValue(const std::string &name,JBinaryDataDef::TpData type){
  const char met[]="CheckSetValue";
  //ResetData();
  int idx=GetValueIndex(name);
  if(idx<0&&GetArray(name)!=NULL)RunException(met,string("The value ")+name+" is an array.");
  if(idx<0&&GetItem(name)!=NULL)RunException(met,string("The value ")+name+" is an item.");
  if(idx>=0&&Values[idx].type!=type)RunException(met,string("Type of value ")+name+" invalid.");
  if(idx<0){
    StValue v;
    if(name.length()>120)RunException(met,string("The name of value ")+name+"  is too large.");
    if(name.empty())RunException(met,string("The name of value ")+name+"  is empty.");
    ResetValue(name,type,v);
    Values.push_back(v);
    idx=int(Values.size())-1;
  }
  ValuesModif=true;
  return(idx);
}

//==============================================================================
/// Reset del valor pasado como referencia.
//==============================================================================
void JBinaryData::ResetValue(const std::string &name,JBinaryDataDef::TpData type,JBinaryData::StValue &v){
  v.name=name; v.type=type; v.vdouble3=TDouble3(0);
}

//==============================================================================
/// Devuelve Value en formato XML.
//==============================================================================
std::string JBinaryData::ValueToXml(const StValue &v)const{
  const char met[]="ValueToXml";
  string tx=JBinaryDataDef::TypeToStr(v.type);
  if(tx.empty())RunException(met,"Name of type invalid.");
  tx=string("<")+tx+" name=\""+v.name+"\" ";
  switch(v.type){
    case JBinaryDataDef::DatText:     tx=tx+"v=\""+ v.vtext +"\" />";                         break;
    case JBinaryDataDef::DatBool:     tx=tx+"v=\""+ (v.vint? "1": "0") +"\" />";              break;    
    case JBinaryDataDef::DatChar:     tx=tx+"v=\""+ fun::IntStr(v.vchar)  +"\" />";           break;    
    case JBinaryDataDef::DatUchar:    tx=tx+"v=\""+ fun::UintStr(v.vuchar)  +"\" />";         break;    
    case JBinaryDataDef::DatShort:    tx=tx+"v=\""+ fun::IntStr(v.vshort) +"\" />";           break;
    case JBinaryDataDef::DatUshort:   tx=tx+"v=\""+ fun::UintStr(v.vushort) +"\" />";         break;
    case JBinaryDataDef::DatInt:      tx=tx+"v=\""+ fun::IntStr(v.vint) +"\" />";             break;
    case JBinaryDataDef::DatUint:     tx=tx+"v=\""+ fun::UintStr(v.vuint) +"\" />";           break;
    case JBinaryDataDef::DatLlong:    tx=tx+"v=\""+ fun::LongStr(v.vllong) +"\" />";          break;
    case JBinaryDataDef::DatUllong:   tx=tx+"v=\""+ fun::UlongStr(v.vullong) +"\" />";        break;
    case JBinaryDataDef::DatFloat:    tx=tx+"v=\""+ fun::FloatStr(v.vfloat,FmtFloat.c_str()) +"\" />";     break;
    case JBinaryDataDef::DatDouble:   tx=tx+"v=\""+ fun::DoubleStr(v.vdouble,FmtDouble.c_str()) +"\" />";  break;
    case JBinaryDataDef::DatInt3:     tx=tx+"x=\""+ fun::IntStr(v.vint3.x) +"\" y=\""+ fun::IntStr(v.vint3.y) +"\" z=\""+ fun::IntStr(v.vint3.z) +"\" />";                                   break;
    case JBinaryDataDef::DatUint3:    tx=tx+"x=\""+ fun::UintStr(v.vuint3.x) +"\" y=\""+ fun::UintStr(v.vuint3.y) +"\" z=\""+ fun::UintStr(v.vuint3.z) +"\" />";                             break;
    case JBinaryDataDef::DatFloat3:   tx=tx+"x=\""+ fun::FloatStr(v.vfloat3.x,FmtFloat.c_str()) +"\" y=\""+ fun::FloatStr(v.vfloat3.y,FmtFloat.c_str()) +"\" z=\""+ fun::FloatStr(v.vfloat3.z,FmtFloat.c_str()) +"\" />";           break;  //"%.7E"
    case JBinaryDataDef::DatDouble3:  tx=tx+"x=\""+ fun::DoubleStr(v.vdouble3.x,FmtDouble.c_str()) +"\" y=\""+ fun::DoubleStr(v.vdouble3.y,FmtDouble.c_str()) +"\" z=\""+ fun::DoubleStr(v.vdouble3.z,FmtDouble.c_str()) +"\" />";  break;  //"%.15E"
    default: RunException(met,"Type of value invalid.");
  }
  return(tx);
}



//==============================================================================
/// Extrae datos del ptr indicado.
//==============================================================================
void JBinaryData::OutData(unsigned &count,unsigned size,const byte *ptr,byte *dat,unsigned sdat)const{
  const unsigned count2=count+sdat;
  if(count2>size)RunException("OutData","Overflow in reading data.");
  memcpy(dat,ptr+count,sdat);
  count=count2;
}

//==============================================================================
/// Extrae string de ptr.
//==============================================================================
std::string JBinaryData::OutStr(unsigned &count,unsigned size,const byte *ptr)const{
  unsigned len=OutUint(count,size,ptr);
  string tex;
  tex.resize(len);
  const unsigned count2=count+len;
  if(count2>size)RunException("OutStr","Overflow in reading data.");
  memcpy((char*)tex.c_str(),ptr+count,len);
  count=count2;
  return(tex);
}

//==============================================================================
/// Introduce datos en ptr.
//==============================================================================
void JBinaryData::InData(unsigned &count,unsigned size,byte *ptr,const byte *dat,unsigned sdat)const{
  const char met[]="InData";
  if(ptr){
    if(count+sdat>size)RunException(met,"Insufficient memory for data.");
    memcpy(ptr+count,dat,sdat);
  }
  if(count+sdat<count)RunException(met,"Size of data is too huge.");
  count+=sdat;
}

//==============================================================================
/// Introduce string en ptr.
//==============================================================================
void JBinaryData::InStr(unsigned &count,unsigned size,byte *ptr,const std::string &cad)const{
  InUint(count,size,ptr,(unsigned)cad.length());
  InData(count,size,ptr,(byte*)cad.c_str(),unsigned(cad.length()));
}

//==============================================================================
/// Introduce Value en ptr.
//==============================================================================
void JBinaryData::InValue(unsigned &count,unsigned size,byte *ptr,const StValue &v)const{
  InStr(count,size,ptr,v.name);
  InInt(count,size,ptr,int(v.type));
  switch(v.type){
    case JBinaryDataDef::DatText:     InStr    (count,size,ptr,v.vtext);     break;
    case JBinaryDataDef::DatBool:     InBool   (count,size,ptr,v.vint!=0);   break;    
    case JBinaryDataDef::DatChar:     InChar   (count,size,ptr,v.vchar);     break;    
    case JBinaryDataDef::DatUchar:    InUchar  (count,size,ptr,v.vuchar);    break;    
    case JBinaryDataDef::DatShort:    InShort  (count,size,ptr,v.vshort);    break;    
    case JBinaryDataDef::DatUshort:   InUshort (count,size,ptr,v.vushort);   break;    
    case JBinaryDataDef::DatInt:      InInt    (count,size,ptr,v.vint);      break;    
    case JBinaryDataDef::DatUint:     InUint   (count,size,ptr,v.vuint);     break;    
    case JBinaryDataDef::DatLlong:    InLlong  (count,size,ptr,v.vllong);    break;     
    case JBinaryDataDef::DatUllong:   InUllong (count,size,ptr,v.vullong);   break;    
    case JBinaryDataDef::DatFloat:    InFloat  (count,size,ptr,v.vfloat);    break;    
    case JBinaryDataDef::DatDouble:   InDouble (count,size,ptr,v.vdouble);   break;    
    case JBinaryDataDef::DatInt3:     InInt3   (count,size,ptr,v.vint3);     break;    
    case JBinaryDataDef::DatUint3:    InUint3  (count,size,ptr,v.vuint3);    break;    
    case JBinaryDataDef::DatFloat3:   InFloat3 (count,size,ptr,v.vfloat3);   break;    
    case JBinaryDataDef::DatDouble3:  InDouble3(count,size,ptr,v.vdouble3);  break;    
    default: RunException("InValue","Type of value invalid.");
  }
}

//==============================================================================
/// Extrae Value en ptr.
//==============================================================================
void JBinaryData::OutValue(unsigned &count,unsigned size,const byte *ptr){
  string name=OutStr(count,size,ptr);
  JBinaryDataDef::TpData type=(JBinaryDataDef::TpData)OutInt(count,size,ptr);
  switch(type){
    case JBinaryDataDef::DatText:     SetvText   (name,OutStr    (count,size,ptr));   break;
    case JBinaryDataDef::DatBool:     SetvBool   (name,OutBool   (count,size,ptr));   break;
    case JBinaryDataDef::DatChar:     SetvChar   (name,OutChar   (count,size,ptr));   break;
    case JBinaryDataDef::DatUchar:    SetvUchar  (name,OutUchar  (count,size,ptr));   break;
    case JBinaryDataDef::DatShort:    SetvShort  (name,OutShort  (count,size,ptr));   break;
    case JBinaryDataDef::DatUshort:   SetvUshort (name,OutUshort (count,size,ptr));   break;
    case JBinaryDataDef::DatInt:      SetvInt    (name,OutInt    (count,size,ptr));   break;
    case JBinaryDataDef::DatUint:     SetvUint   (name,OutUint   (count,size,ptr));   break;
    case JBinaryDataDef::DatLlong:    SetvLlong  (name,OutLlong  (count,size,ptr));   break;
    case JBinaryDataDef::DatUllong:   SetvUllong (name,OutUllong (count,size,ptr));   break;
    case JBinaryDataDef::DatFloat:    SetvFloat  (name,OutFloat  (count,size,ptr));   break;
    case JBinaryDataDef::DatDouble:   SetvDouble (name,OutDouble (count,size,ptr));   break;
    case JBinaryDataDef::DatInt3:     SetvInt3   (name,OutInt3   (count,size,ptr));   break;
    case JBinaryDataDef::DatUint3:    SetvUint3  (name,OutUint3  (count,size,ptr));   break;
    case JBinaryDataDef::DatFloat3:   SetvFloat3 (name,OutFloat3 (count,size,ptr));   break;
    case JBinaryDataDef::DatDouble3:  SetvDouble3(name,OutDouble3(count,size,ptr));   break;
    default: RunException("OutValue","Type of value invalid.");
  }
}

//==============================================================================
/// Introduce datos basicos de Array en ptr.
//==============================================================================
void JBinaryData::InArrayBase(unsigned &count,unsigned size,byte *ptr,const JBinaryDataArray *ar)const{
  InStr(count,size,ptr,CodeArrayDef);
  InStr(count,size,ptr,ar->GetName());
  InBool(count,size,ptr,ar->GetHide());
  InInt(count,size,ptr,int(ar->GetType()));
  InUint(count,size,ptr,ar->GetCount());
  //-Calcula e introduce size de los datos del array.
  unsigned sizearraydata=0;
  InArrayData(sizearraydata,0,NULL,ar);
  InUint(count,size,ptr,sizearraydata);
}
//==============================================================================
/// Introduce contendido de Array en ptr.
//==============================================================================
void JBinaryData::InArrayData(unsigned &count,unsigned size,byte *ptr,const JBinaryDataArray *ar)const{
  const JBinaryDataDef::TpData type=ar->GetType();
  const unsigned num=ar->GetCount();
  const void* pointer=ar->GetPointer();
  if(num&&!pointer)RunException("InArrayData","Pointer of array with data is invalid.");
  if(type==JBinaryDataDef::DatText){//-Array de strings.
    const string *list=(string*)pointer;
    for(unsigned c=0;c<num;c++)InStr(count,size,ptr,list[c]);
  }
  else{//-Array de tipos basicos.
    unsigned sizetype=(unsigned)JBinaryDataDef::SizeOfType(ar->GetType());
    InData(count,size,ptr,(byte*)pointer,sizetype*num);
  }
}

//==============================================================================
/// Introduce Array en ptr.
//==============================================================================
void JBinaryData::InArray(unsigned &count,unsigned size,byte *ptr,const JBinaryDataArray *ar)const{
  //-Calcula size de la definicion del array.
  unsigned sizearray=0;
  InArrayBase(sizearray,0,NULL,ar);
  //-Introduce propiedades de array en ptr.
  InUint(count,size,ptr,sizearray);
  InArrayBase(count,size,ptr,ar);
  //-Introduce contenido del array.
  InArrayData(count,size,ptr,ar);
}

//==============================================================================
/// Introduce datos basicos de Item en ptr.
//==============================================================================
void JBinaryData::InItemBase(unsigned &count,unsigned size,byte *ptr,bool all)const{
  InStr(count,size,ptr,CodeItemDef);
  InStr(count,size,ptr,GetName());
  InBool(count,size,ptr,GetHide());     
  InBool(count,size,ptr,GetHideValues());
  InStr(count,size,ptr,GetFmtFloat());
  InStr(count,size,ptr,GetFmtDouble());
  InUint(count,size,ptr,(all? GetArraysCount(): GetVisibleArraysCount()));
  InUint(count,size,ptr,(all? GetItemsCount(): GetVisibleItemsCount()));
  if(all||!HideValues)InUint(count,size,ptr,(ValuesModif? GetSizeValues(): ValuesSize));
  else InUint(count,size,ptr,0);
}

//==============================================================================
/// Introduce Item en ptr.
//==============================================================================
void JBinaryData::InItem(unsigned &count,unsigned size,byte *ptr,bool all)const{
  //-Calcula size de la definicion del item.
  unsigned sizeitem=0;
  InItemBase(sizeitem,0,NULL,all);
  //-Introduce propiedades de item en ptr.
  InUint(count,size,ptr,sizeitem);
  InItemBase(count,size,ptr,all);
  //-Introduce values en ptr.
  if(all||!HideValues){
    if(ValuesModif)SaveValues(count,size,ptr);//-Cache no valida.
    else InData(count,size,ptr,ValuesData,ValuesSize);//-Cache actualizada.
  }
  //-Introduce arrays en ptr.
  for(unsigned c=0;c<Arrays.size();c++)if(all||!Arrays[c]->GetHide())InArray(count,size,ptr,Arrays[c]);
  //-Introduce items en ptr.
  for(unsigned c=0;c<Items.size();c++)if(all||!Items[c]->GetHide())Items[c]->InItem(count,size,ptr,all);
}

//==============================================================================
/// Extrae datos basicos del Array de ptr.
//==============================================================================
JBinaryDataArray* JBinaryData::OutArrayBase(unsigned &count,unsigned size,const byte *ptr,unsigned &countdata,unsigned &sizedata){
  const char met[]="OutArrayBase";
  if(OutStr(count,size,ptr)!=CodeArrayDef)RunException(met,"Validation code is invalid.");
  string name=OutStr(count,size,ptr);
  bool hide=OutBool(count,size,ptr);
  JBinaryDataDef::TpData type=(JBinaryDataDef::TpData)OutInt(count,size,ptr);
  countdata=OutUint(count,size,ptr);
  sizedata=OutUint(count,size,ptr);
  if(type!=JBinaryDataDef::DatText&&sizedata!=JBinaryDataDef::SizeOfType(type)*countdata)RunException(met,"Size of data is invalid.");
  //-Crea array.
  JBinaryDataArray *ar=CreateArray(name,type);
  ar->SetHide(hide);
  return(ar);
}

//==============================================================================
/// Extrae contenido de Array de ptr.
//==============================================================================
void JBinaryData::OutArrayData(unsigned &count,unsigned size,const byte *ptr,JBinaryDataArray *ar,unsigned countdata,unsigned sizedata){
  if(ar->GetType()==JBinaryDataDef::DatText){//-Array de strings.
    ar->AllocMemory(countdata);
    for(unsigned c=0;c<countdata;c++)ar->AddText(OutStr(count,size,ptr),false);
  }
  else{
    //-Comprueba que los datos del array estan disponibles.
    unsigned count2=count+sizedata;
    if(count2>size)RunException("OutArrayData","Overflow in reading data.");
    //-Extrae datos para el array.
    ar->AddData(countdata,ptr+count,true);
    count=count2;
  }
}

//==============================================================================
/// Extrae Array de ptr.
//==============================================================================
void JBinaryData::OutArray(unsigned &count,unsigned size,const byte *ptr){
  //-Crea y configura array a partir de ptr.
  const unsigned sizearraydef=OutUint(count,size,ptr);
  unsigned countdata,sizedata;
  JBinaryDataArray *ar=OutArrayBase(count,size,ptr,countdata,sizedata);
  //-Extrae contenido del array.
  OutArrayData(count,size,ptr,ar,countdata,sizedata);
}

//==============================================================================
/// Extrae propiedades basicas de Item de ptr.
//==============================================================================
JBinaryData* JBinaryData::OutItemBase(unsigned &count,unsigned size,const byte *ptr,bool create,unsigned &narrays,unsigned &nitems,unsigned &sizevalues){
  if(OutStr(count,size,ptr)!=CodeItemDef)RunException("OutItem","Validation code is invalid.");
  JBinaryData* item=this;
  if(create)item=CreateItem(OutStr(count,size,ptr));
  else item->SetName(OutStr(count,size,ptr));
  item->SetHide(OutBool(count,size,ptr));
  item->SetHideValues(OutBool(count,size,ptr),false);
  item->SetFmtFloat(OutStr(count,size,ptr),false);
  item->SetFmtDouble(OutStr(count,size,ptr),false);
  narrays=OutUint(count,size,ptr);
  nitems=OutUint(count,size,ptr);
  sizevalues=OutUint(count,size,ptr);
  return(item);
}

//==============================================================================
/// Extrae Item de ptr.
//==============================================================================
void JBinaryData::OutItem(unsigned &count,unsigned size,const byte *ptr,bool create){
  //-Extrae propiedades del item.
  const unsigned sizeitemdef=OutUint(count,size,ptr);
  unsigned narrays,nitems,sizevalues;
  JBinaryData* item=OutItemBase(count,size,ptr,create,narrays,nitems,sizevalues);
  //-Extrae values del item.
  if(sizevalues){
    if(OutStr(count,size,ptr)!=CodeValuesDef)RunException("OutItem","Validation code is invalid.");
    unsigned num=OutUint(count,size,ptr);
    for(unsigned c=0;c<num;c++)item->OutValue(count,size,ptr);
  }
  //-Extrae arrays del item.
  for(unsigned c=0;c<narrays;c++)item->OutArray(count,size,ptr);
  //-Extrae items del item.
  for(unsigned c=0;c<nitems;c++)item->OutItem(count,size,ptr,true);
}

//==============================================================================
/// Devuelve el tama�o necesario para almacenar todos los values del item.
//==============================================================================
unsigned JBinaryData::GetSizeValues()const{
  unsigned count=0;
  SaveValues(count,0,NULL);
  return(count);
}

//==============================================================================
/// Almacena datos de values en ptr.
//==============================================================================
void JBinaryData::SaveValues(unsigned &count,unsigned size,byte *ptr)const{
  unsigned num=unsigned(Values.size());
  InStr(count,size,ptr,CodeValuesDef);
  InUint(count,size,ptr,num);
  for(unsigned c=0;c<num;c++)InValue(count,size,ptr,Values[c]);
}



//==============================================================================
/// Graba contenido de array en fichero.
//==============================================================================
void JBinaryData::WriteArrayData(std::fstream *pf,const JBinaryDataArray *ar)const{
  const JBinaryDataDef::TpData type=ar->GetType();
  const unsigned countdata=ar->GetCount();
  const void* pointer=ar->GetPointer();
  if(countdata&&!pointer)RunException("WriteArrayData","Pointer of array with data is invalid.");
  if(type==JBinaryDataDef::DatText){//-Array de strings.
    const string *list=(string*)pointer;
    unsigned sbuf=0;
    for(unsigned c=0;c<countdata;c++)InStr(sbuf,0,NULL,list[c]);//-Calcula size de buffer.
    byte *buf=new byte[sbuf];
    unsigned cbuf=0;
    for(unsigned c=0;c<countdata;c++)InStr(cbuf,sbuf,buf,list[c]);//-Copia en buffer.
    pf->write((char*)buf,cbuf);
    delete[] buf;
  }
  else{//-Array de tipos basicos.
    unsigned sizetype=(unsigned)JBinaryDataDef::SizeOfType(ar->GetType());
    pf->write((char*)pointer,sizetype*countdata);
  }
}

//==============================================================================
/// Graba Array en fichero.
//==============================================================================
void JBinaryData::WriteArray(std::fstream *pf,unsigned sbuf,byte *buf,const JBinaryDataArray *ar)const{
  //-Calcula size de la definicion del array.
  unsigned sizearray=0;
  InArrayBase(sizearray,0,NULL,ar);
  //-Graba propiedades de array.
  unsigned cbuf=0;
  InUint(cbuf,sbuf,buf,sizearray);
  InArrayBase(cbuf,sbuf,buf,ar);
  pf->write((char*)buf,cbuf);
  //-Graba contenido del array.
  WriteArrayData(pf,ar);
}

//==============================================================================
/// Graba Item en fichero.
//==============================================================================
void JBinaryData::WriteItem(std::fstream *pf,unsigned sbuf,byte *buf,bool all)const{
  //-Calcula size de la definicion del item.
  unsigned sizeitem=0;
  InItemBase(sizeitem,0,NULL,all);
  //-Graba propiedades de item en fichero.
  {
    unsigned cbuf=0;
    InUint(cbuf,sbuf,buf,sizeitem);
    InItemBase(cbuf,sbuf,buf,all);
    pf->write((char*)buf,cbuf);
  }
  //-Graba values.
  if(all||!GetHideValues())pf->write((char*)ValuesData,ValuesSize);
  //-Graba arrays.
  for(unsigned c=0;c<Arrays.size();c++)if(all||!Arrays[c]->GetHide())WriteArray(pf,sbuf,buf,Arrays[c]);
  //-Graba items.
  for(unsigned c=0;c<Items.size();c++)if(all||!Items[c]->GetHide())Items[c]->WriteItem(pf,sbuf,buf,all);
}


//==============================================================================
/// Devuelve unsigned leido de fichero.
//==============================================================================
unsigned JBinaryData::ReadUint(std::ifstream *pf)const{
  unsigned v=0;
  pf->read((char*)&v,sizeof(unsigned));
  return(v);
}

//==============================================================================
/// Carga datos de array de fichero.
//==============================================================================
void JBinaryData::ReadArrayData(std::ifstream *pf,JBinaryDataArray *ar,unsigned countdata,unsigned sizedata,bool loadarraysdata){
  const JBinaryDataDef::TpData type=ar->GetType();
  if(loadarraysdata)ar->ReadData(countdata,sizedata,pf,true);
  else{
    ar->ConfigFileData((llong)pf->tellg(),countdata,sizedata);  
    pf->seekg(sizedata,ios::cur);
  }
}

//==============================================================================
/// Carga array de fichero.
//==============================================================================
void JBinaryData::ReadArray(std::ifstream *pf,unsigned sbuf,byte *buf,bool loadarraysdata){
  //-Carga propiedades del array.
  const unsigned sizearraydef=ReadUint(pf);
  pf->read((char*)buf,sizearraydef);
  unsigned countdata,sizedata;
  unsigned cbuf=0;
  JBinaryDataArray *ar=OutArrayBase(cbuf,sizearraydef,buf,countdata,sizedata);
  //-Extrae contenido del array.
  ReadArrayData(pf,ar,countdata,sizedata,loadarraysdata);
}

//==============================================================================
/// Carga Item de fichero.
//==============================================================================
void JBinaryData::ReadItem(std::ifstream *pf,unsigned sbuf,byte *buf,bool create,bool loadarraysdata){
  //-Carga propiedades del item.
  const unsigned sizeitemdef=ReadUint(pf);
  pf->read((char*)buf,sizeitemdef);
  unsigned narrays,nitems,sizevalues;
  unsigned cbuf=0;
  JBinaryData* item=OutItemBase(cbuf,sizeitemdef,buf,create,narrays,nitems,sizevalues);
  //-Carga values del item.
  if(sizevalues){
    byte* buf2=(sbuf>=sizevalues? buf: NULL);
    if(!buf2)buf2=new byte[sizevalues];
    pf->read((char*)buf2,sizevalues);
    unsigned cbuf2=0;
    if(OutStr(cbuf2,sizevalues,buf2)!=CodeValuesDef)RunException("ReadItem","Validation code is invalid.");
    unsigned num=OutUint(cbuf2,sizevalues,buf2);
    for(unsigned c=0;c<num;c++)item->OutValue(cbuf2,sizevalues,buf2);
    if(buf!=buf2)delete[] buf2;
  }
  //-Carga arrays del item.
  for(unsigned c=0;c<narrays;c++)item->ReadArray(pf,sbuf,buf,loadarraysdata);
  //-Carga items del item.
  for(unsigned c=0;c<nitems;c++)item->ReadItem(pf,sbuf,buf,true,loadarraysdata);
}

//==============================================================================
/// Genera cabecera para fichero.
//==============================================================================
JBinaryData::StHeadFmtBin JBinaryData::MakeFileHead(const std::string &filecode)const{
  StHeadFmtBin hfmt; 
  memset(&hfmt,0,sizeof(StHeadFmtBin));
  string titu=string("#FileJBD ")+filecode;
  unsigned stitu=min(58u,unsigned(titu.size()));
  for(unsigned c=0;c<stitu;c++)hfmt.titu[c]=titu[c];
  for(unsigned c=stitu;c<58;c++)hfmt.titu[c]=' ';
  hfmt.titu[58]='\n';
  hfmt.byteorder=byte(fun::GetByteOrder());
  return(hfmt);
}

//==============================================================================
/// Devuelve tama�o de fichero y su cabecera.
/// Si el fichero no contiene una cabecera devuelve 0.
//==============================================================================
unsigned JBinaryData::GetFileHead(std::ifstream *pf,JBinaryData::StHeadFmtBin &head)const{
  //-Obtiene size del fichero.
  pf->seekg(0,ios::end);
  const unsigned fsize=(unsigned)pf->tellg();   //printf("CheckFileHead> FileSize:%u\n",fsize);
  pf->seekg(0,ios::beg);
  //-Lee cabecera basica.
  if(fsize>=sizeof(StHeadFmtBin))pf->read((char*)&head,sizeof(StHeadFmtBin));
  else memset(&head,0,sizeof(StHeadFmtBin));
  return(fsize);
}

//==============================================================================
/// Comprueba formato de cabecera con filecode y bitorder.
/// Genera excepcion en caso de error.
//==============================================================================
void JBinaryData::CheckHead(const std::string &file,const StHeadFmtBin &head,const std::string &filecode)const{
  const char met[]="CheckHead";
  int err=0;
  //-Coprueba formato de cabecera y filecode.
  if(!err){
    StHeadFmtBin head2=MakeFileHead(filecode);
    unsigned c=0;
    for(;head.titu[c]==head2.titu[c]&&c<60;c++);
    if(c<9)err=2;
    else if(!filecode.empty()&&c<60)err=3;
  }
  //-Coprueba orden de bytes.
  if(!err){
    byte byteorder=byte(fun::GetByteOrder());
    if(head.byteorder!=byte(fun::GetByteOrder()))err=1;
  }
  if(err==1)RunException(met,"The byte-order in file is invalid.",file);
  else if(err==2)RunException(met,"The format file JBinaryData is invalid.",file);
  else if(err==3)RunException(met,"The file code is invalid.",file);
}

//==============================================================================
/// Comprueba formato de cabecera con filecode y bitorder.
/// Si el fichero este vacio tambien genera excepcion.
/// Devuelve size de fichero.
//==============================================================================
unsigned JBinaryData::CheckFileHead(const std::string &file,std::ifstream *pf,const std::string &filecode)const{
  JBinaryData::StHeadFmtBin head;
  //-Obtiene size y cabecera del fichero.
  const unsigned fsize=GetFileHead(pf,head);
  //-Comprueba validez de cabecera.
  CheckHead(file,head,filecode);
  return(fsize);
}

//==============================================================================
/// Comprueba formato de cabecera con filecode y bitorder.
/// En caso de que el fichero este vacio no genera excepcion.
/// Devuelve size de fichero.
//==============================================================================
unsigned JBinaryData::CheckFileListHead(const std::string &file,std::fstream *pf,const std::string &filecode)const{
  const char met[]="CheckFileHead";
  //-Obtiene size del fichero.
  pf->seekg(0,ios::end);
  const unsigned fsize=(unsigned)pf->tellg();   //printf("CheckFileHead> FileSize:%u\n",fsize);
  pf->seekg(0,ios::beg);
  //-Lee cabecera basica y comprueba validez.
  StHeadFmtBin head;
  if(fsize>=sizeof(StHeadFmtBin)){
    pf->read((char*)&head,sizeof(StHeadFmtBin));
    //-Comprueba validez de cabecera.
    CheckHead(file,head,filecode);
  }
  return(fsize);
}

//==============================================================================
/// Graba contenido en fichero XML.
//==============================================================================
void JBinaryData::WriteFileXmlArray(const std::string &tabs,std::ofstream* pf,bool svarrays,const JBinaryDataArray* ar)const{
  const char met[]="WriteFileXmlArray";
  const JBinaryDataDef::TpData type=ar->GetType();
  const unsigned size=ar->GetSize();
  const unsigned count=ar->GetCount();
  string tx=JBinaryDataDef::TypeToStr(type);
  if(tx.empty())RunException(met,"Name of type invalid.");
  string res=string("<array_")+tx+" name=\""+ar->GetName()+"\" size=\""+fun::UintStr(size)+"\" count=\""+fun::UintStr(count)+"\" hide=\""+ (ar->GetHide()? '1': '0') +"\"";
  if(!svarrays){
    res=res+"/>";
    (*pf) << tabs << res << endl;
  }
  else{
    res=res+">";
    (*pf) << tabs << res << endl;
    const void *data=ar->GetDataPointer();
    JBinaryData::StValue v;
    //ResetValue("",type,v);
    v.type=type;
    for(unsigned c=0;c<count;c++){
      v.name=fun::UintStr(c);
      switch(v.type){
        case JBinaryDataDef::DatText:     v.vtext=   ((const string *)       data)[c];   break;
        case JBinaryDataDef::DatBool:     v.vint=   (((const bool *)data)[c]? 1: 0);     break;
        case JBinaryDataDef::DatChar:     v.vchar=   ((const char *)         data)[c];   break;
        case JBinaryDataDef::DatUchar:    v.vuchar=  ((const unsigned char *)data)[c];   break;
        case JBinaryDataDef::DatShort:    v.vshort=  ((const short *)        data)[c];   break;
        case JBinaryDataDef::DatUshort:   v.vushort= ((const word *)         data)[c];   break;
        case JBinaryDataDef::DatInt:      v.vint=    ((const int *)          data)[c];   break;   
        case JBinaryDataDef::DatUint:     v.vuint=   ((const unsigned *)     data)[c];   break;   
        case JBinaryDataDef::DatLlong:    v.vllong=  ((const llong *)        data)[c];   break;   
        case JBinaryDataDef::DatUllong:   v.vullong= ((const ullong *)       data)[c];   break;   
        case JBinaryDataDef::DatFloat:    v.vfloat=  ((const float *)        data)[c];   break;   
        case JBinaryDataDef::DatDouble:   v.vdouble= ((const double *)       data)[c];   break;   
        case JBinaryDataDef::DatInt3:     v.vint3=   ((const tint3 *)        data)[c];   break;   
        case JBinaryDataDef::DatUint3:    v.vuint3=  ((const tuint3 *)       data)[c];   break;   
        case JBinaryDataDef::DatFloat3:   v.vfloat3= ((const tfloat3 *)      data)[c];   break;   
        case JBinaryDataDef::DatDouble3:  v.vdouble3=((const tdouble3 *)     data)[c];   break;   
        default: RunException(met,"Type of value invalid.");
      }
      (*pf) << tabs << "\t" << ValueToXml(v) << endl;
    }
    (*pf) << tabs << string("</array_")+tx+">" << endl;
  }
}

//==============================================================================
/// Graba contenido en fichero XML.
//==============================================================================
void JBinaryData::WriteFileXml(const std::string &tabs,std::ofstream* pf,bool svarrays)const{
  (*pf) << tabs << "<item name=\"" << GetName() << "\" hide=\"" << (GetHide()? '1': '0') << "\" hidevalues=\"" << (GetHideValues()? '1': '0') << "\">" << endl;
  for(unsigned c=0;c<Values.size();c++)(*pf) << tabs << "\t" << ValueToXml(Values[c]) << endl;
  for(unsigned c=0;c<Arrays.size();c++)WriteFileXmlArray(tabs+"\t",pf,svarrays,Arrays[c]);
  for(unsigned c=0;c<Items.size();c++)Items[c]->WriteFileXml(tabs+"\t",pf,svarrays);
  (*pf) << tabs << "</item>" << endl;
}

//==============================================================================
/// Cambia nombre de objeto comprobando que no tenga hermanos con el mismo nombre.
//==============================================================================
void JBinaryData::SetName(const std::string &name){
  const char met[]="SetName";
  if(Parent){
    if(Parent->ExistsValue(name))RunException(met,"There is already a value with the name given.");
    if(Parent->GetArray(name)!=NULL)RunException(met,"There is already an array with the name given.");
    if(Parent->GetItem(name)!=NULL)RunException(met,"There is already an item with the name given.");
  }
  Name=name;
}

//==============================================================================
/// Cambia oculatacion de values.
//==============================================================================
void JBinaryData::SetHideValues(bool hide,bool down){
  HideValues=hide;
  if(down)for(unsigned c=0;c<Items.size();c++)Items[c]->SetHideValues(hide,true);
}

//==============================================================================
/// Cambia oculatacion de arrays.
//==============================================================================
void JBinaryData::SetHideArrays(bool hide,bool down){
  for(unsigned c=0;c<Arrays.size();c++){
    Arrays[c]->SetHide(hide);
    if(down)for(unsigned c=0;c<Items.size();c++)Items[c]->SetHideArrays(hide,true);
  }
}

//==============================================================================
/// Cambia oculatacion de items.
//==============================================================================
void JBinaryData::SetHideItems(bool hide,bool down){
  for(unsigned c=0;c<Items.size();c++){
    Items[c]->SetHide(hide);
    if(down)Items[c]->SetHideItems(hide,true);
  }
}

//==============================================================================
/// Cambia formato de texto para valores float.
//==============================================================================
void JBinaryData::SetFmtFloat(const std::string &fmt,bool down){
  FmtFloat=fmt;
  if(down)for(unsigned c=0;c<Items.size();c++)Items[c]->SetFmtFloat(fmt,true);
}

//==============================================================================
/// Cambia formato de texto para valores double.
//==============================================================================
void JBinaryData::SetFmtDouble(const std::string &fmt,bool down){
  FmtDouble=fmt;
  if(down)for(unsigned c=0;c<Items.size();c++)Items[c]->SetFmtDouble(fmt,true);
}

//==============================================================================
/// Devuelve el tama�o necesario para almacenar todos los datos del item y descendientes.
/// Con all activado se incluyen tambien los elementos ocultos.
//==============================================================================
unsigned JBinaryData::GetSizeDataConst(bool all)const{
  unsigned count=0;
  InItem(count,0,NULL,all);
  return(count);
}

//==============================================================================
/// Almacena datos de item en ptr y devuelve los bytes almacenados.
/// Con all activado se incluyen tambien los elementos ocultos.
//==============================================================================
unsigned JBinaryData::SaveDataConst(unsigned size,byte* ptr,bool all)const{
  if(!ptr)RunException("SaveDataConst","The pointer is invalid.");
  unsigned count=0;
  InItem(count,size,ptr,all);
  return(count);
}

//==============================================================================
/// Devuelve el tama�o necesario para almacenar todos los del item y 
/// descendientes. Actualiza la cache de values  para mejorar el rendimiento en 
/// operaciones posteriores.
/// Con all activado se incluyen tambien los elementos ocultos.
//==============================================================================
unsigned JBinaryData::GetSizeData(bool all){
  ValuesCachePrepare(true);
  unsigned count=0;
  InItem(count,0,NULL,all);
  return(count);
}

//==============================================================================
/// Almacena datos de item en ptr y devuelve los bytes almacenados. Actualiza la
/// cache de values para mejorar el rendimiento en operaciones posteriores.
/// Con all activado se incluyen tambien los elementos ocultos.
//==============================================================================
unsigned JBinaryData::SaveData(unsigned size,byte* ptr,bool all){
  if(!ptr)RunException("SaveData","The pointer is invalid.");
  ValuesCachePrepare(true);
  unsigned count=0;
  InItem(count,size,ptr,all);
  return(count);
}

//==============================================================================
/// Carga datos desde ptr.
//==============================================================================
void JBinaryData::LoadData(unsigned size,const byte* ptr){
  if(!ptr)RunException("LoadData","The pointer is invalid.");
  Clear(); //-Limpia contenido de objeto.
  unsigned count=0;
  OutItem(count,size,ptr,false);
}


//==============================================================================
/// Graba datos en fichero.
/// Con memory utiliza un buffer para todos los datos. Consume mas memoria pero
/// es mas rapido si no hay arrays grandes.
/// Con all activado se incluyen tambien los elementos ocultos.
//==============================================================================
void JBinaryData::SaveFileData(std::fstream *pf,bool head,const std::string &filecode,bool memory,bool all)const{
  if(head){//-Graba cabecera basica.
    StHeadFmtBin head=MakeFileHead(filecode); 
    pf->write((char*)&head,sizeof(StHeadFmtBin));
  }
  //-Graba datos.
  if(memory){//-Graba datos desde memoria.
    const unsigned sbuf=GetSizeDataConst(all);
    //printf("SaveFile> sbuf:%u\n",sbuf);
    byte *buf=new byte[sbuf];
    unsigned sbuf2=SaveDataConst(sbuf,buf,all);
    //printf("SaveFile> sbuf2:%u\n",sbuf2);
    pf->write((char*)buf,sbuf);
    delete[] buf;
  }
  else{//-Graba datos directamente.
    const unsigned sbuf=1024;
    byte buf[sbuf];
    WriteItem(pf,sbuf,buf,all);
  }
}

//==============================================================================
/// Graba datos en fichero.
/// Con memory utiliza un buffer para todos los datos. Consume mas memoria pero
/// es mas rapido si no hay arrays grandes.
/// Con all activado se incluyen tambien los elementos ocultos.
//==============================================================================
void JBinaryData::SaveFile(const std::string &file,bool memory,bool all){
  const char met[]="SaveFile";
  ValuesCachePrepare(true);
  fstream pf;
  pf.open(file.c_str(),ios::binary|ios::out);
  if(pf){
    SaveFileData(&pf,true,Name,memory,all);
    if(pf.fail())RunException(met,"File writing failure.",file);
    pf.close();
  }
  else RunException(met,"Cannot open the file.",file);
}

//==============================================================================
/// Carga datos de un fichero.
/// Con memory utiliza un buffer para todos los datos. Consume mas memoria pero
/// es mas rapido cuando no hay arrays grandes.
//==============================================================================
void JBinaryData::LoadFile(const std::string &file,const std::string &filecode,bool memory){
  const char met[]="LoadFile";
  Clear(); //-Limpia contenido de objeto.
  ifstream pf;
  pf.open(file.c_str(),ios::binary|ios::in);
  if(pf){
    const unsigned fsize=CheckFileHead(file,&pf,filecode);
    //-Carga datos.
    if(memory){//-Carga datos desde memoria.
      const unsigned sbuf=fsize-sizeof(StHeadFmtBin);
      //printf("LoadFile> sbuf:%u\n",sbuf);
      byte *buf=new byte[sbuf];
      pf.read((char*)buf,sbuf);
      LoadData(sbuf,buf);
      delete[] buf;
    }
    else{//-Carga datos directamente.
      const unsigned sbuf=1024;
      byte buf[sbuf];
      ReadItem(&pf,sbuf,buf,false,true);
    }
    pf.close();
  }
  else RunException(met,"Cannot open the file.",file);
}

//==============================================================================
/// Graba Item en fichero a continuacion de los items que ya exitan.
/// Con memory utiliza un buffer para todos los datos. Consume mas memoria pero
/// es mas rapido si no hay arrays grandes.
/// Con all activado se incluyen tambien los elementos ocultos.
//==============================================================================
void JBinaryData::SaveFileListApp(const std::string &file,const std::string &filecode,bool memory,bool all){
  const char met[]="SaveFile";
  ValuesCachePrepare(true);
  fstream pf;
  if(fun::FileExists(file))pf.open(file.c_str(),ios::binary|ios::out|ios::in|ios::app);
  else pf.open(file.c_str(),ios::binary|ios::out);
  if(pf){
    const unsigned fsize=CheckFileListHead(file,&pf,filecode);
    pf.seekp(0,pf.end);
    //-Graba datos de parent.
    if(!fsize)Parent->SaveFileData(&pf,true,filecode,memory,all);
    //-Graba datos de item.
    SaveFileData(&pf,false,filecode,memory,all);
    if(pf.fail())RunException(met,"File writing failure.",file);
    pf.close();
  }
  else RunException(met,"Cannot open the file.",file);
}

//==============================================================================
/// Carga datos de un fichero.
/// Con memory utiliza un buffer para todos los datos. Consume mas memoria pero
/// es mas rapido cuando no hay arrays grandes.
//==============================================================================
void JBinaryData::LoadFileListApp(const std::string &file,const std::string &filecode,bool memory){
  const char met[]="LoadFileListApp";
  Clear(); //-Limpia contenido de objeto.
  ifstream pf;
  pf.open(file.c_str(),ios::binary|ios::in);
  if(pf){
    SetName(filecode);
    const unsigned fsize=CheckFileHead(file,&pf,filecode);
    unsigned pfile=(unsigned)pf.tellg();
    while(pfile<fsize){
      //-Carga datos.
      if(memory){//-Carga datos desde memoria.
        const unsigned sbuf=fsize-sizeof(StHeadFmtBin);
        //printf("LoadFile> sbuf:%u\n",sbuf);
        byte *buf=new byte[sbuf];
        pf.read((char*)buf,sbuf);
        unsigned cbuf=0;
        while(cbuf<sbuf){
          OutItem(cbuf,sbuf,buf,true);
          //-Renombra ultimo item leido.
          unsigned lastitem=GetItemsCount()-1;
          JBinaryData* ite=GetItem(lastitem);
          ite->SetName(fun::PrintStr("LS%04u_",lastitem)+ite->GetName());
        }
        delete[] buf;
      }
      else{//-Carga datos directamente.
        const unsigned sbuf=1024;
        byte buf[sbuf];
        ReadItem(&pf,sbuf,buf,true,true);
        //-Renombra ultimo item leido.
        unsigned lastitem=GetItemsCount()-1;
        JBinaryData* ite=GetItem(lastitem);
        ite->SetName(fun::PrintStr("LS%04u_",lastitem)+ite->GetName());
      }
      pfile=(unsigned)pf.tellg();
    }
    pf.close();
  }
  else RunException(met,"Cannot open the file.",file);
}

//==============================================================================
/// Abre fichero y carga estructura de datos pero sin cargar el contenido de los
/// arrays.
//==============================================================================
void JBinaryData::OpenFileStructure(const std::string &file,const std::string &filecode){
  const char met[]="OpenFileStructure";
  if(Parent)RunException(met,"Item is not root.");
  Clear(); //-Limpia contenido de objeto.
  FileStructure=new ifstream;
  FileStructure->open(file.c_str(),ios::binary|ios::in);
  if(*FileStructure){
    const unsigned fsize=CheckFileHead(file,FileStructure,filecode);
    const unsigned sbuf=1024;
    byte buf[sbuf];
    ReadItem(FileStructure,sbuf,buf,false,false);
  }
  else{
    CloseFileStructure();
    RunException(met,"Cannot open the file.",file);
  }
}

//==============================================================================
/// Cierra fichero abierto con OpenFileStructure().
//==============================================================================
void JBinaryData::CloseFileStructure(){
  if(FileStructure&&FileStructure->is_open())FileStructure->close();
  delete FileStructure; FileStructure=NULL;
}

//==============================================================================
/// Devuelve puntero al fichero abierto con OpenFileStructure().
//==============================================================================
std::ifstream* JBinaryData::GetFileStructure()const{
  if(Parent)RunException("GetFileStructure","Item is not root.");
  return(FileStructure);
}


//==============================================================================
/// Graba contenido en fichero XML.
//==============================================================================
void JBinaryData::SaveFileXml(std::string file,bool svarrays,const std::string &head)const{
  const char* met="SaveFileXml";
  file=fun::GetWithoutExtension(file)+".xml";
  ofstream pf;
  pf.open(file.c_str());
  if(pf){
    pf << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << endl;
    pf << "<data" << head <<" date=\"" << fun::GetDateTime() << "\">" << endl;
    WriteFileXml("\t",&pf,svarrays);
    pf << "</data>" << endl;
    if(pf.fail())RunException(met,"Failed writing to file.",file);
    pf.close();
  }
  else RunException(met,"File could not be opened.",file);
}




//==============================================================================
/// Devuelve item principal.
//==============================================================================
JBinaryData* JBinaryData::GetItemRoot(){
  return(Parent? Parent->GetItemRoot(): this);
}

//==============================================================================
/// Devuelve el numero de items no marcados como ocultos.
//==============================================================================
unsigned JBinaryData::GetVisibleItemsCount()const{
  unsigned num=0;
  for(unsigned c=0;c<Items.size();c++)if(!Items[c]->GetHide())num++;
  return(num);
}

//==============================================================================
/// Devuelve indice del item con el nombre indicado o -1 si no existe.
//==============================================================================
int JBinaryData::GetItemIndex(const std::string &name){
  int idx=-1;
  for(unsigned c=0;c<Items.size()&&idx<0;c++)if(Items[c]->Name==name)idx=c;
  return(idx);
}

//==============================================================================
/// Devuelve item con el nombre indicado o NULL si no existe.
//==============================================================================
JBinaryData* JBinaryData::GetItem(const std::string &name){
  JBinaryData* ret=NULL;
  for(unsigned c=0;c<Items.size()&&!ret;c++)if(Items[c]->Name==name)ret=Items[c];
  return(ret);
}

//==============================================================================
/// Devuelve item segun el indice indicado o NULL si no existe.
//==============================================================================
JBinaryData* JBinaryData::GetItem(unsigned index){
  return(index>=GetItemsCount()? NULL: Items[index]);
}

//==============================================================================
/// Crea y devuelve item con el nombre. Genera excepcion si ya existe.
//==============================================================================
JBinaryData* JBinaryData::CreateItem(const std::string &name){
  if(GetItem(name)!=NULL)RunException("CreateItem","There is already an item with the name given.");
  JBinaryData *item=new JBinaryData(name);
  item->Parent=this;
  Items.push_back(item);
  return(item);
}

//==============================================================================
/// Elimina el item indicado.
//==============================================================================
void JBinaryData::RemoveItem(const std::string &name){
  int idx=GetItemIndex(name);
  if(idx>=0){
    JBinaryData* item=Items[idx];
    Items.erase(Items.begin()+idx);
    delete item;
  }
}

//==============================================================================
/// Elimina todos los items almacenados.
//==============================================================================
void JBinaryData::RemoveItems(){
  for(unsigned c=0;c<Items.size();c++)delete Items[c];
  Items.clear();
}




//==============================================================================
/// Devuelve el numero de arrays no marcados como ocultos.
//==============================================================================
unsigned JBinaryData::GetVisibleArraysCount()const{
  unsigned num=0;
  for(unsigned c=0;c<Arrays.size();c++)if(!Arrays[c]->GetHide())num++;
  return(num);
}

//==============================================================================
/// Devuelve posicion de la variable solicitada, -1 en caso de no existir.
//==============================================================================
int JBinaryData::GetArrayIndex(const std::string &name)const{
  int idx=-1; 
  for(unsigned c=0;c<Arrays.size()&&idx<0;c++)if(Arrays[c]->GetName()==name)idx=c;
  return(idx);
}

//==============================================================================
/// Devuelve array con el nombre indicado o NULL si no existe.
//==============================================================================
JBinaryDataArray* JBinaryData::GetArray(const std::string &name){
  JBinaryDataArray* ret=NULL;
  for(unsigned c=0;c<Arrays.size()&&!ret;c++)if(Arrays[c]->GetName()==name)ret=Arrays[c];
  return(ret);
}

//==============================================================================
/// Devuelve array segun el indice indicado o NULL si no existe.
//==============================================================================
JBinaryDataArray* JBinaryData::GetArray(unsigned index){
  return(index>=GetArraysCount()? NULL: Arrays[index]);
}

//==============================================================================
/// Crea y devuelve array con el nombre. Genera excepcion si ya existe.
//==============================================================================
JBinaryDataArray* JBinaryData::CreateArray(const std::string &name,JBinaryDataDef::TpData type){
  if(GetItem(name)!=NULL)RunException("CreateArray","There is already an array with the name given.");
  JBinaryDataArray *ar=new JBinaryDataArray(this,name,type);
  Arrays.push_back(ar);
  return(ar);
}

//==============================================================================
/// Crea y devuelve array con datos.
//==============================================================================
JBinaryDataArray* JBinaryData::CreateArray(const std::string &name,JBinaryDataDef::TpData type,unsigned count,const void *data,bool externalpointer){
  JBinaryDataArray *ar=CreateArray(name,type);
  ar->SetData(count,data,externalpointer);
  return(ar);
}

//==============================================================================
/// Elimina el array indicado.
//==============================================================================
void JBinaryData::RemoveArray(const std::string &name){
  int idx=GetArrayIndex(name);
  if(idx>=0){
    JBinaryDataArray* ar=Arrays[idx];
    Arrays.erase(Arrays.begin()+idx);
    delete ar;
  }
}

//==============================================================================
/// Elimina todos los arrays almacenados.
//==============================================================================
void JBinaryData::RemoveArrays(){
  for(unsigned c=0;c<Arrays.size();c++)delete Arrays[c];
  Arrays.clear();
}




//==============================================================================
/// Devuelve posicion de la variable solicitada, -1 en caso de no existir.
//==============================================================================
int JBinaryData::GetValueIndex(const std::string &name)const{
  int pos=-1; 
  for(unsigned c=0;c<Values.size()&&pos<0;c++)if(Values[c].name==name)pos=c;
  return(pos);
}

//==============================================================================
/// Devuelve el nombre del value solicitado (vacio si no existe).
//==============================================================================
std::string JBinaryData::NameOfValue(unsigned index)const{
  return(index>=GetValuesCount()? "": Values[index].name);
}

//==============================================================================
/// Devuelve el tipo del value solicitado (DatNull si no existe).
//==============================================================================
JBinaryDataDef::TpData JBinaryData::TypeOfValue(const std::string &name)const{
  int idx=GetValueIndex(name);
  return(idx<0? JBinaryDataDef::DatNull: Values[idx].type);
}

//==============================================================================
/// Devuelve el tipo del value solicitado (DatNull si no existe).
//==============================================================================
JBinaryDataDef::TpData JBinaryData::TypeOfValue(unsigned index)const{
  return(index>=GetValuesCount()? JBinaryDataDef::DatNull: Values[index].type);
}

//==============================================================================
/// Indica si existe el value solicitado.
//==============================================================================
bool JBinaryData::ExistsValue(const std::string &name)const{
  return(GetValueIndex(name)>=0);
}

//==============================================================================
/// Indica si existe el value solicitado del tipo indicado.
//==============================================================================
bool JBinaryData::ExistsValue(const std::string &name,JBinaryDataDef::TpData type)const{
  int idx=GetValueIndex(name);
  return(idx>=0&&Values[idx].type==type);
}

//==============================================================================
/// Elimina el value indicado.
//==============================================================================
void JBinaryData::RemoveValue(const std::string &name){
  int idx=GetValueIndex(name);
  if(idx>=0)Values.erase(Values.begin()+idx);
  ValuesModif=true;
}

//==============================================================================
/// Elimina todos los values almacenados.
//==============================================================================
void JBinaryData::RemoveValues(){
  Values.clear();
  ValuesCacheReset();
}



//==============================================================================
/// Devuelve el valor solicitado de tipo texto.
//==============================================================================
std::string JBinaryData::GetvText(const std::string &name,bool optional,std::string valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatText);
  return(pos<0? valdef: Values[pos].vtext);
}
//==============================================================================
/// Devuelve el valor solicitado de tipo bool.
//==============================================================================
bool JBinaryData::GetvBool(const std::string &name,bool optional,bool valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatBool);
  return(pos<0? valdef: Values[pos].vint!=0);
}
//==============================================================================
/// Devuelve el valor solicitado de tipo char.
//==============================================================================
char JBinaryData::GetvChar(const std::string &name,bool optional,char valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatChar);
  return(pos<0? valdef: Values[pos].vchar);
}
//==============================================================================
/// Devuelve el valor solicitado de tipo unsigned char.
//==============================================================================
unsigned char JBinaryData::GetvUchar(const std::string &name,bool optional,unsigned char valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatUchar);
  return(pos<0? valdef: Values[pos].vuchar);
}
//==============================================================================
/// Devuelve el valor solicitado de tipo short.
//==============================================================================
short JBinaryData::GetvShort(const std::string &name,bool optional,short valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatShort);
  return(pos<0? valdef: Values[pos].vshort);
}
//==============================================================================
/// Devuelve el valor solicitado de tipo unsigned short.
//==============================================================================
unsigned short JBinaryData::GetvUshort(const std::string &name,bool optional,unsigned short valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatUshort);
  return(pos<0? valdef: Values[pos].vushort);
}
//==============================================================================
/// Devuelve el valor solicitado de tipo int.
//==============================================================================
int JBinaryData::GetvInt(const std::string &name,bool optional,int valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatInt);
  return(pos<0? valdef: Values[pos].vint);
}
//==============================================================================
/// Devuelve el valor solicitado de tipo unsigned.
//==============================================================================
unsigned JBinaryData::GetvUint(const std::string &name,bool optional,unsigned valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatUint);
  return(pos<0? valdef: Values[pos].vuint);
}
//==============================================================================
/// Devuelve el valor solicitado de tipo long long.
//==============================================================================
llong JBinaryData::GetvLlong(const std::string &name,bool optional,llong valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatLlong);
  return(pos<0? valdef: Values[pos].vllong);
}
//==============================================================================
/// Devuelve el valor solicitado de tipo unsigned long long.
//==============================================================================
ullong JBinaryData::GetvUllong(const std::string &name,bool optional,ullong valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatUllong);
  return(pos<0? valdef: Values[pos].vullong);
}
//==============================================================================
/// Devuelve el valor solicitado de tipo float.
//==============================================================================
float JBinaryData::GetvFloat(const std::string &name,bool optional,float valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatFloat);
  return(pos<0? valdef: Values[pos].vfloat);
}
//==============================================================================
/// Devuelve el valor solicitado de tipo double.
//==============================================================================
double JBinaryData::GetvDouble(const std::string &name,bool optional,double valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatDouble);
  return(pos<0? valdef: Values[pos].vdouble);
}
//==============================================================================
/// Devuelve el valor solicitado de tipo tint3.
//==============================================================================
tint3 JBinaryData::GetvInt3(const std::string &name,bool optional,tint3 valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatInt3);
  return(pos<0? valdef: Values[pos].vint3);
}
//==============================================================================
/// Devuelve el valor solicitado de tipo tuint3.
//==============================================================================
tuint3 JBinaryData::GetvUint3(const std::string &name,bool optional,tuint3 valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatUint3);
  return(pos<0? valdef: Values[pos].vuint3);
}
//==============================================================================
/// Devuelve el valor solicitado de tipo tfloat3.
//==============================================================================
tfloat3 JBinaryData::GetvFloat3(const std::string &name,bool optional,tfloat3 valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatFloat3);
  return(pos<0? valdef: Values[pos].vfloat3);
}
//==============================================================================
/// Devuelve el valor solicitado de tipo tdouble3.
//==============================================================================
tdouble3 JBinaryData::GetvDouble3(const std::string &name,bool optional,tdouble3 valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatDouble3);
  return(pos<0? valdef: Values[pos].vdouble3);
}



//==============================================================================
/// Crea o modifica un valor de tipo texto.
//==============================================================================
void JBinaryData::SetvText(const std::string &name,const std::string &v){
  Values[CheckSetValue(name,JBinaryDataDef::DatText)].vtext=v;
}
//==============================================================================
/// Crea o modifica un valor de tipo bool.
//==============================================================================
void JBinaryData::SetvBool(const std::string &name,bool v){
  Values[CheckSetValue(name,JBinaryDataDef::DatBool)].vint=(v? 1: 0);
}
//==============================================================================
/// Crea o modifica un valor de tipo char.
//==============================================================================
void JBinaryData::SetvChar(const std::string &name,char v){
  Values[CheckSetValue(name,JBinaryDataDef::DatChar)].vchar=v;
}
//==============================================================================
/// Crea o modifica un valor de tipo unsigned char.
//==============================================================================
void JBinaryData::SetvUchar(const std::string &name,unsigned char v){
  Values[CheckSetValue(name,JBinaryDataDef::DatUchar)].vuchar=v;
}
//==============================================================================
/// Crea o modifica un valor de tipo short.
//==============================================================================
void JBinaryData::SetvShort(const std::string &name,short v){
  Values[CheckSetValue(name,JBinaryDataDef::DatShort)].vshort=v;
}
//==============================================================================
/// Crea o modifica un valor de tipo unsigned short.
//==============================================================================
void JBinaryData::SetvUshort(const std::string &name,unsigned short v){
  Values[CheckSetValue(name,JBinaryDataDef::DatUshort)].vushort=v;
}
//==============================================================================
/// Crea o modifica un valor de tipo int.
//==============================================================================
void JBinaryData::SetvInt(const std::string &name,int v){
  Values[CheckSetValue(name,JBinaryDataDef::DatInt)].vint=v;
}
//==============================================================================
/// Crea o modifica un valor de tipo unsigned.
//==============================================================================
void JBinaryData::SetvUint(const std::string &name,unsigned v){
  Values[CheckSetValue(name,JBinaryDataDef::DatUint)].vuint=v;
}
//==============================================================================
/// Crea o modifica un valor de tipo long long.
//==============================================================================
void JBinaryData::SetvLlong(const std::string &name,llong v){
  Values[CheckSetValue(name,JBinaryDataDef::DatLlong)].vllong=v;
}
//==============================================================================
/// Crea o modifica un valor de tipo unsigned long long.
//==============================================================================
void JBinaryData::SetvUllong(const std::string &name,ullong v){
  Values[CheckSetValue(name,JBinaryDataDef::DatUllong)].vullong=v;
}
//==============================================================================
/// Crea o modifica un valor de tipo float.
//==============================================================================
void JBinaryData::SetvFloat(const std::string &name,float v){
  Values[CheckSetValue(name,JBinaryDataDef::DatFloat)].vfloat=v;
}
//==============================================================================
/// Crea o modifica un valor de tipo double.
//==============================================================================
void JBinaryData::SetvDouble(const std::string &name,double v){
  Values[CheckSetValue(name,JBinaryDataDef::DatDouble)].vdouble=v;
}
//==============================================================================
/// Crea o modifica un valor de tipo tint3.
//==============================================================================
void JBinaryData::SetvInt3(const std::string &name,tint3 v){
  Values[CheckSetValue(name,JBinaryDataDef::DatInt3)].vint3=v;
}
//==============================================================================
/// Crea o modifica un valor de tipo tuint3.
//==============================================================================
void JBinaryData::SetvUint3(const std::string &name,tuint3 v){
  Values[CheckSetValue(name,JBinaryDataDef::DatUint3)].vuint3=v;
}
//==============================================================================
/// Crea o modifica un valor de tipo tfloat3.
//==============================================================================
void JBinaryData::SetvFloat3(const std::string &name,tfloat3 v){
  Values[CheckSetValue(name,JBinaryDataDef::DatFloat3)].vfloat3=v;
}
//==============================================================================
/// Crea o modifica un valor de tipo tdouble3.
//==============================================================================
void JBinaryData::SetvDouble3(const std::string &name,tdouble3 v){
  Values[CheckSetValue(name,JBinaryDataDef::DatDouble3)].vdouble3=v;
}




