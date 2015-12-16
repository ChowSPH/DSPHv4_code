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

#include "JReadDatafile.h"
#include "Functions.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JReadDatafile::JReadDatafile(){
  ClassName="JReadDatafile";
  Data=NULL; LineBegin=NULL;
  Reset();
}
//==============================================================================
/// Destructor.
//==============================================================================
JReadDatafile::~JReadDatafile(){
  Reset();
}
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JReadDatafile::Reset(){
  delete[] Data; Data=NULL;
  delete[] LineBegin; LineBegin=NULL;
  SizeFile=Size=0;
  File=""; Sep=";";
  LineCount=RemLineCount=0;
  ResetReadLine();
}

//==============================================================================
/// Load file data.
//==============================================================================
void JReadDatafile::LoadFile(const std::string &file,unsigned maxsize){
  const char met[]="LoadFile";
  Reset();
  File=file;
  ifstream pf;
  pf.open(file.c_str(),ios::binary);
  if(pf){
    pf.seekg(0,ios::end);
    SizeFile=unsigned(pf.tellg())+1;
    if(SizeFile>maxsize)RunException(met,"File exceeds the maximum size allowed.",file);
    Data=new char[SizeFile];
    pf.seekg(0,ios::beg);
    pf.read(Data,SizeFile-1);
    Data[SizeFile-1]='\n';
    pf.close();
  }
  else RunException(met,"Cannot open the file.",file);
  Size=SizeFile;
  ProcessLines();
  SetReadLine(0);
}

//==============================================================================
/// Replaces several spaces and tabulations by a tabulation.
/// Removes spaces and tabulations at begin and end of line.
//==============================================================================
void JReadDatafile::ProcessSpaces(){
  unsigned nsp=0;
  unsigned c2=0;
  bool begin=true;
  for(unsigned c=0;c<Size;c++){
    char let=Data[c];
    if(let==' ' || let=='\t')nsp++;
    else{
      if(let=='\n'){ nsp=0; begin=true; }
      if(nsp){ 
        if(!begin){ Data[c2]='\t'; c2++; }
        nsp=0; 
        begin=false;
      }
      if(c!=c2)Data[c2]=let;
      c2++;
    }
  }
  //if(nsp){ Data[c2]='\t'; c2++; nsp=0; }
  //printf("++> Remove spaces+tabs: %u -> %u\n",Size,c2);
  Size=c2;
}

//==============================================================================
/// Load begin and end of each line and counts lines.
//==============================================================================
void JReadDatafile::ProcessLines(){
  const char met[]="MakeLineBegin";
  //-Remove character \r and counts lines.
  LineCount=0;
  unsigned c2=0;
  for(unsigned c=0;c<Size;c++){
    char let=Data[c];
    if(let=='\n')LineCount++;
    if(c!=c2)Data[c2]=let;
    if(let!='\r')c2++;
  }
  //printf("++> Remove \\r: %u -> %u\n",Size,c2);
  Size=c2;
  //-Remove empty lines at end.
  while(Size>1 && Data[Size-2]=='\n'){ Size--; LineCount--; }
  //printf("++> Remove empty lines: %u\n",c2-Size);
  //-Allocate memory.
  LineBegin=new unsigned[LineCount+1];
  //-Prepares lines and looks for separator.
  bool run=true;
  while(run){
    run=false;
    //-Load LineBegin[].
    unsigned lin=0;
    LineBegin[0]=0;
    for(unsigned c=0;c<Size;c++){
      char let=Data[c];
      if(let=='\n'){
        lin++; LineBegin[lin]=c+1;
      }
    }
    if(lin!=LineCount)RunException(met,"Error counting lines.");
    LineBegin[LineCount]=Size;
    //-Counts remark lines.
    RemLineCount=0;
    for(int c=0;c<LineCount;c++)if(Data[LineBegin[c]]=='#')RemLineCount++;
    //printf("++> RemLineCount: %u\n",RemLineCount);
    //-Prints lines.
    if(0)for(int c=0;c<LineCount;c++){
      const unsigned pini=LineBegin[c];
      const unsigned pfin=LineBegin[c+1];
      if(!c)printf("\n");
      printf("++> Line[%02d]=[%2d]=[",c,pfin-pini-1); 
      for(unsigned p=pini;p<pfin-1;p++)printf("%c",Data[p]);
      printf("]\n");
    }
    //-Determines the separator.
    {
      unsigned sep0=0,sep1=0,sep2=0,sep3=0;
      unsigned nlin=20;
      for(int c=0;c<LineCount && nlin;c++)if(Data[LineBegin[c]]!='#'){
        nlin--;
        const unsigned pini=LineBegin[c];
        const unsigned pfin=LineBegin[c+1];
        for(unsigned p=pini;p<pfin;p++){
          switch(Data[p]){
            case ' ':   sep0++;  break;
            case '\t':  sep1++;  break;
            case ';':   sep2++;  break;
            case ',':   sep3++;  break;
          }
        }
      }
      //printf("++> Sep: [ ]:%u [\\t]:%u [,]:%u [;]:%u\n",sep0,sep1,sep2,sep3);
      sep1+=sep0;
      if(sep1>=sep2 && sep1>=sep3)Sep="\t";
      else if(sep2>=sep1 && sep2>=sep3)Sep=";";
      else Sep=",";
      if(Sep=="\t" && sep0){
        ProcessSpaces();
        run=true;
      }
    }
  }
  //throw "finito...";
  //printf("++> LineCount:%u(%u)  Size:%u -> %u\n",LineCount,RemLineCount,SizeFile,Size);
}

//==============================================================================
/// Restarts reading position.
//==============================================================================
void JReadDatafile::ResetReadLine(){ 
  ReadLin=ReadLinValue=-1;
  ReadLine="";
  ReadValue="";
}

//==============================================================================
/// Configures reading position. (line=0...n-1)
//==============================================================================
void JReadDatafile::SetReadLine(int line){ 
  ResetReadLine();
  ReadLin=line;
  ReadLine=GetLine(ReadLin);
  ReadLinValue=-1;
}

//==============================================================================
/// Returns selected line.
//==============================================================================
std::string JReadDatafile::GetLine(int line)const{
  if(line<0 || line>=LineCount)RunException("GetLine","Line number is invalid.",File);
  const unsigned ini=LineBegin[line];
  const unsigned len=LineBegin[line+1]-ini-1;//-Con -1 se quita el salto de linea.
  string tex;
  if(len){
    tex.resize(len);
    memcpy((char*)tex.c_str(),Data+ini,len);
  }
  return(tex);
}

//==============================================================================
/// Returns next value in the current line or next line if in_line is false.
//==============================================================================
std::string JReadDatafile::ReadNextValue(bool in_line){
  if(in_line && ReadLine.empty())RunException("ReadNextValue",fun::PrintStr("Value %d does not exist in line %d.",ReadLinValue+1,ReadLin),File);
  while(ReadLine.empty() || ReadLine[0]=='#')SetReadLine(ReadLin+1);
  //printf("==>> ReadLine:[%s] Sep:[%u]\n",ReadLine.c_str(),Sep[0]);
  ReadValue=fun::StrSplit(Sep,ReadLine); 
  ReadLinValue++;
  return(ReadValue);
}

//==============================================================================
/// Returns next double in the current line or next line if in_line is false.
//==============================================================================
double JReadDatafile::ReadNextDouble(bool in_line){
  string value=ReadNextValue(in_line);
  return(atof(ReadValue.c_str()));
}

//==============================================================================
/// Returns next tdouble3 in the current line or next line if in_line is false.
//==============================================================================
tdouble3 JReadDatafile::ReadNextDouble3(bool in_line){
  tdouble3 v;
  v.x=ReadNextDouble(in_line);
  v.y=ReadNextDouble(true);
  v.z=ReadNextDouble(true);
  return(v);
}

//==============================================================================
/// Returns next int in the current line or next line if in_line is false.
//==============================================================================
int JReadDatafile::ReadNextInt(bool in_line){
  string value=ReadNextValue(in_line);
  return(atoi(ReadValue.c_str()));
}


