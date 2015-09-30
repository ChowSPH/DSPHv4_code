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

#include "JRadixSort.h"
#include <string>
#include <cstring>
#include <climits>
using namespace std;

#ifdef _WITHOMP
  #include <omp.h>  //Activar tb en Props config -> C/C++ -> Lenguaje -> OpenMp
#else
  #define omp_get_thread_num() 0
  #define omp_get_max_threads() 1
#endif


//==============================================================================
/// Constructor de objetos.
//==============================================================================
JRadixSort::JRadixSort(bool useomp):UseOmp(useomp){
  ClassName="JRadixSort";
  InitData32=NULL; InitData64=NULL;
  Data32=NULL; Data64=NULL;
  PrevData32=NULL; PrevData64=NULL;
  BeginKeys=NULL;
  Index=NULL; PrevIndex=NULL;
  Reset();
}

//==============================================================================
/// Destructor de objetos.
//==============================================================================
JRadixSort::~JRadixSort(){
  Reset();
}

//==============================================================================
/// Reinicializa el estado del objeto, recuperando la configuracion por defecto.
//==============================================================================
void JRadixSort::Reset(){
  if(Data32==InitData32)Data32=NULL;
  if(Data64==InitData64)Data64=NULL;
  if(PrevData32==InitData32)PrevData32=NULL;
  if(PrevData64==InitData64)PrevData64=NULL;
  InitData32=NULL; InitData64=NULL;
  delete[] Data32; Data32=NULL;
  delete[] Data64; Data64=NULL;
  delete[] PrevData32; PrevData32=NULL;
  delete[] PrevData64; PrevData64=NULL;
  Size=Nbits=Nkeys=0;
  delete[] BeginKeys; BeginKeys=NULL;
  delete[] Index; Index=NULL;
  delete[] PrevIndex; PrevIndex=NULL;
}

//==============================================================================
/// Devuelve el numero de bits necesarios para codificar el valor indicado.
//==============================================================================
template<class T> unsigned JRadixSort::TBitsSize(T v,unsigned smax)const{
  unsigned nbits=1;
  for(;v>>nbits&&nbits<smax;nbits++);
/*
  printf("++++> v>>nbits(%u): %u\n",nbits,(v>>nbits));
  for(;v>>nbits&&nbits<smax;nbits++){
    printf("++++> v>>nbits(%u): %u\n",nbits,(v>>nbits));
    //if(nbits>65)break;
  }
  printf("++++> BitSize(%u):%u\n",v,nbits);
*/
  return(nbits);
}
//==============================================================================
unsigned JRadixSort::BitsSize(unsigned v)const{ return(TBitsSize<unsigned>(v,32)); };
unsigned JRadixSort::BitsSize(unsigned long long v)const{ return(TBitsSize<unsigned long long>(v,64)); };

//==============================================================================
/// Calcula Nbits para los datos indicados.
//==============================================================================
template<class T> unsigned JRadixSort::TCalcNbits(unsigned size,const T *data)const{
  const char met[]="CalcNbits";
  const int threads=omp_get_max_threads();
  T mxdata=0;

  if(!UseOmp || threads<2){//-Secuencial.
    T vmax=0;
    for(unsigned c=0;c<size;c++)vmax=max(vmax,data[c]);
    mxdata=vmax;
  }
  else{//-Con OpenMP.
    //-Calcula bloques de ejecucion.
    const int nk=int(size/OMPSIZE)+1;
    if(nk<0)RunException(met,"Number of values is invalid.");
    const int rk=int(size%OMPSIZE);
    //printf("CalcNbits> Size:%u  nk:%d  rk:%d  sizek:%u\n",size,nk,rk,OMPSIZE);
    //-Calcula maximo de nk bloques con varios hilos.
    T *vmax=new T[threads*OMPSTRIDE];
    memset(vmax,0,sizeof(T)*threads*OMPSTRIDE);
    #ifdef _WITHOMP
      #pragma omp parallel for schedule (static)
    #endif
    for(int c=0;c<nk;c++){
      int th=omp_get_thread_num();
      T mx=vmax[OMPSTRIDE*th];
      const unsigned c2ini=OMPSIZE*c;
      const unsigned c2fin=c2ini+(c+1<nk? OMPSIZE: rk);
      for(unsigned c2=c2ini;c2<c2fin;c2++)mx=max(mx,data[c2]);
      vmax[OMPSTRIDE*th]=mx;
    }
    //-Calcula reduce maximo de todos los hilos.
    T mx=0;
    for(int t=0;t<threads;t++)mx=max(mx,vmax[OMPSTRIDE*t]);
    delete[] vmax; vmax=NULL;
    mxdata=mx;
  }
  //-Calcula nbits para el valor maximo.
  return(BitsSize(mxdata));
}
//==============================================================================
unsigned JRadixSort::CalcNbits(unsigned size,const unsigned *data)const{ return(TCalcNbits<unsigned>(size,data)); }
unsigned JRadixSort::CalcNbits(unsigned size,const unsigned long long *data)const{ return(TCalcNbits<unsigned long long>(size,data)); }

//==============================================================================
/// Reserva la memoria necesaria.
//==============================================================================
void JRadixSort::AllocMemory(unsigned s){
  try{
    if(Type32)Data32=new unsigned[s];
    else Data64=new unsigned long long[s];
  }
  catch(const std::bad_alloc){
    RunException("AllocMemory","Cannot allocate the requested memory.");
  }
}

//==============================================================================
/// Contabiliza numero de valores para cada clave.
//==============================================================================
template<class T> void JRadixSort::LoadBeginKeys(const T* data){
  const char met[]="LoadBeginKeys";
  const int threads=omp_get_max_threads();
  //-Reserva espacio para contadores de claves.
  Nkeys=unsigned((Nbits+(KEYSBITS-1))/KEYSBITS);
  delete[] BeginKeys; BeginKeys=NULL;
  BeginKeys=new unsigned[Nkeys*KEYSRANGE];

  //-Inicia proceso.
  if(!UseOmp || threads<2){//-Secuencial.
    unsigned *nkeys=new unsigned[Nkeys*KEYSRANGE];
    memset(nkeys,0,sizeof(unsigned)*Nkeys*KEYSRANGE);
    for(unsigned c2=0;c2<Size;c2++){
      const T v=data[c2];
      for(unsigned ck=0;ck<Nkeys;ck++){ 
        const unsigned k=((v>>(ck*KEYSBITS))&KEYSMASK);
        nkeys[ck*KEYSRANGE+k]++;
      } 
    }
    //-Carga valores en BeginKeys.
    for(unsigned ck=0;ck<Nkeys;ck++){
      BeginKeys[ck*KEYSRANGE]=0;
      for(unsigned c=1;c<KEYSRANGE;c++){
        BeginKeys[ck*KEYSRANGE+c]=BeginKeys[ck*KEYSRANGE+c-1]+nkeys[ck*KEYSRANGE+c-1];
        //printf(">> BeginKeys[%u*%u+%u]=%u\n",ck,KEYSRANGE,c,BeginKeys[ck*KEYSRANGE+c]);
      }
    }
    delete[] nkeys;
  }
  else{//-con OpenMP.
    //-Calcula bloques de ejecucion.
    const int nk=int(Size/OMPSIZE)+1;
    if(nk<0)RunException(met,"Number of values is invalid.");
    const int rk=int(Size%OMPSIZE);
    //-Reserva memoria auxiliar para conteo.
    const unsigned skeys=Nkeys*KEYSRANGE+100;
    unsigned *nkeys=new unsigned[skeys*threads];
    memset(nkeys,0,sizeof(unsigned)*skeys*threads);
    //-Realiza conteo con varios hilos.
    //printf(">> nk:%d  rk:%d\n",nk,rk);
    #ifdef _WITHOMP
      #pragma omp parallel for schedule (static)
    #endif
    for(int c=0;c<nk;c++){
      int th=omp_get_thread_num();
      unsigned *n=nkeys+(skeys*th);
      const unsigned c2ini=OMPSIZE*c;
      const unsigned c2fin=c2ini+(c+1<nk? OMPSIZE: rk);
      //printf(">> c2ini:%d  c2fin:%d\n",c2ini,c2fin);
      for(unsigned c2=c2ini;c2<c2fin;c2++){
        T v=data[c2];
        for(unsigned ck=0;ck<Nkeys;ck++){ 
          unsigned k=((v>>(ck*KEYSBITS))&KEYSMASK);
          n[ck*KEYSRANGE+k]++;
          //printf(">> acc n[%u*KEYSRANGE+%u]=%u\n",ck,c,n[ck*KEYSRANGE+k]);
        } 
      }
    }
    //-Reduce conteo de todos los hilos.
    //for(int t=0;t<Threads;t++)for(unsigned ck=0;ck<Nkeys;ck++)for(unsigned c=0;c<KEYSRANGE;c++)
    //  printf(">> nkeys[%u*skeys+%u*%u+%u]=%u\n",t,ck,KEYSRANGE,c,nkeys[t*skeys+ck*KEYSRANGE+c]);
    for(int t=1;t<threads;t++)for(unsigned ck=0;ck<Nkeys;ck++)for(unsigned c=0;c<KEYSRANGE;c++)nkeys[ck*KEYSRANGE+c]+=nkeys[t*skeys+ck*KEYSRANGE+c];
    //-Carga valores en BeginKeys.
    for(unsigned ck=0;ck<Nkeys;ck++){
      BeginKeys[ck*KEYSRANGE]=0;
      for(unsigned c=1;c<KEYSRANGE;c++){
        BeginKeys[ck*KEYSRANGE+c]=BeginKeys[ck*KEYSRANGE+c-1]+nkeys[ck*KEYSRANGE+c-1];
        //printf(">> BeginKeys[%u*%u+%u]=%u\n",ck,KEYSRANGE,c,BeginKeys[ck*KEYSRANGE+c]);
      }
    }
    //-Libera memoria auxiliar.
    delete[] nkeys;
  }
}

//==============================================================================
/// Realiza un paso de ordenacion en funcion de 1 bit.
//==============================================================================
template<class T> void JRadixSort::SortStep(unsigned ck,const T* data,T* data2){
  const char met[]="SortStep";
  unsigned p2[KEYSRANGE];
  memcpy(p2,BeginKeys+(ck*KEYSRANGE),sizeof(unsigned)*KEYSRANGE);
  const unsigned ckmov=ck*KEYSBITS;
  for(unsigned p=0;p<Size;p++){
    unsigned k=((data[p]>>ckmov)&KEYSMASK);
    data2[p2[k]]=data[p];
    p2[k]++;
  }
  if(0){//dbg
    //printf("CheckSortStep2[%d]>\n",ck);
    //for(unsigned c=0;c<KEYSRANGE-1;c++){
    //  printf("  > p2[%d]: %u != %u :BeginKeys[ck*KEYSRANGE+c+1]\n",c,p2[c],BeginKeys[ck*KEYSRANGE+c+1]);
    //}
    //printf("  > p2[%d]: %u != %u :Size\n",KEYSRANGE-1,p2[KEYSRANGE-1],Size);
    for(unsigned c=0;c<KEYSRANGE-1;c++)if(p2[c]!=BeginKeys[ck*KEYSRANGE+c+1])RunException(met,"El numero de puntos ordenados no es correcto");
    if(p2[KEYSRANGE-1]!=Size)RunException(met,"El numero de puntos ordenados no es correcto");
  }
}

//==============================================================================
/// Realiza un paso de ordenacion en funcion de 1 bit.
//==============================================================================
template<class T> void JRadixSort::SortStepIndex(unsigned ck,const T* data,T* data2,const unsigned *index,unsigned *index2){
  const char met[]="SortStep";
  unsigned p2[KEYSRANGE];
  memcpy(p2,BeginKeys+(ck*KEYSRANGE),sizeof(unsigned)*KEYSRANGE);
  const unsigned ckmov=ck*KEYSBITS;
  for(unsigned p=0;p<Size;p++){
    unsigned k=((data[p]>>ckmov)&KEYSMASK);
    unsigned pk=p2[k];
    data2[pk]=data[p];
    index2[pk]=index[p];
    p2[k]++;
  }
  if(0){//dbg
    //printf("CheckSortStep2[%d]>\n",ck);
    //for(unsigned c=0;c<KEYSRANGE-1;c++){
    //  printf("  > p2[%d]: %u != %u :BeginKeys[ck*KEYSRANGE+c+1]\n",c,p2[c],BeginKeys[ck*KEYSRANGE+c+1]);
    //}
    //printf("  > p2[%d]: %u != %u :Size\n",KEYSRANGE-1,p2[KEYSRANGE-1],Size);
    for(unsigned c=0;c<KEYSRANGE-1;c++)if(p2[c]!=BeginKeys[ck*KEYSRANGE+c+1])RunException(met,"El numero de puntos ordenados no es correcto");
    if(p2[KEYSRANGE-1]!=Size)RunException(met,"El numero de puntos ordenados no es correcto");
  }
}

//==============================================================================
/// Crea e inicializa el vector Index[].
//==============================================================================
void JRadixSort::IndexCreate(){
  const char met[]="IndexCreate";
  const int threads=omp_get_max_threads();
  //-Reserva memoria.
  try{
    Index=new unsigned[Size];
    PrevIndex=new unsigned[Size];
    if(0){ memset(Index,255,sizeof(unsigned)*Size); }//dbg
  }
  catch(const std::bad_alloc){
    RunException(met,"Cannot allocate the requested memory.");
  }

  //-Carga PrevIndex[] con valores consecutivos.
  if(!UseOmp || threads<2){//-Secuencial.
    for(unsigned c2=0;c2<Size;c2++)PrevIndex[c2]=c2;
  }
  else{//-con OpenMP.
    const int nk=int(Size/OMPSIZE)+1;
    if(nk<0)RunException(met,"Number of values is invalid.");
    const int rk=int(Size%OMPSIZE);
    //-Realiza proceso con varios hilos.
    //printf(">> nk:%d  rk:%d\n",nk,rk);
    #ifdef _WITHOMP
      #pragma omp parallel for schedule (static)
    #endif
    for(int c=0;c<nk;c++){
      const unsigned c2ini=OMPSIZE*c;
      const unsigned c2fin=c2ini+(c+1<nk? OMPSIZE: rk);
      //printf(">> c2ini:%d  c2fin:%d\n",c2ini,c2fin);
      for(unsigned c2=c2ini;c2<c2fin;c2++)PrevIndex[c2]=c2;
    }
  }
  //if(1){
  //  bool ok=true;
  //  for(unsigned p=0;p<Size&&ok;p++)ok=(PrevIndex[p]==p);
  //  printf("---IndexCreate>\n");
  //  if(!ok)throw "Error";
  //}
}

//==============================================================================
/// Ordena valores de data.
//==============================================================================
void JRadixSort::Sort(bool makeindex,unsigned size,unsigned *data,unsigned nbits){
  Reset();
  Nbits=nbits; Size=size; 
  Type32=true; InitData32=data; PrevData32=data;
  AllocMemory(Size);
  if(makeindex)IndexCreate();
  LoadBeginKeys<unsigned>(PrevData32);
  //if(1)for(unsigned b=0;b<Nbits;b++)printf("Sort> CountZeros[%d]:%d + %d\n",b,CountZeros[b],Size-CountZeros[b]); //dbg
  for(unsigned ck=0;ck<Nkeys;ck++){
    if(makeindex){
      SortStepIndex(ck,PrevData32,Data32,PrevIndex,Index);
      swap(PrevIndex,Index);
    }
    else SortStep(ck,PrevData32,Data32);
    swap(PrevData32,Data32);
  }
  if(makeindex){ 
    swap(PrevIndex,Index);
    delete[] PrevIndex; PrevIndex=NULL;
  }
  //-Copia los datos en el puntero recibido como parametro.
  if(PrevData32!=InitData32)memcpy(InitData32,PrevData32,sizeof(unsigned)*Size);
}

//==============================================================================
/// Ordena valores de data.
//==============================================================================
void JRadixSort::Sort(bool makeindex,unsigned size,unsigned long long *data,unsigned nbits){
  Reset();
  Nbits=nbits; Size=size; 
  Type32=false; InitData64=data; PrevData64=data;
  AllocMemory(Size);
  if(makeindex)IndexCreate();
  LoadBeginKeys<unsigned long long>(PrevData64);
  //if(1)for(unsigned b=0;b<Nbits;b++)printf("Sort> CountZeros[%d]:%d + %d\n",b,CountZeros[b],Size-CountZeros[b]); //dbg
  for(unsigned ck=0;ck<Nkeys;ck++){
    if(makeindex){
      SortStepIndex(ck,PrevData64,Data64,PrevIndex,Index);
      swap(PrevIndex,Index);
    }
    else SortStep(ck,PrevData64,Data64);
    swap(PrevData64,Data64);
    //printf("\n***** SortStep: %d\n",ck);
    //for(unsigned p=0;p<Size&&p<100;p++){ 
    //  const unsigned ckmov=ck*KEYSBITS;
    //  unsigned k=((PrevData64[p]>>ckmov)&KEYSMASK);
    //  printf("data[%2d]: %12llu  k:%u\n",p,PrevData64[p],k);
    //}
  }
  if(makeindex){ 
    swap(PrevIndex,Index);
    delete[] PrevIndex; PrevIndex=NULL;
  }
  //-Copia los datos en el puntero recibido como parametro.
  if(PrevData64!=InitData64)memcpy(InitData64,PrevData64,sizeof(unsigned long long)*Size);
}

//==============================================================================
/// Crea indice para ordenacion pero sin modificar los datos pasados.
//==============================================================================
void JRadixSort::MakeIndex(unsigned size,const unsigned *data,unsigned nbits){
  unsigned* auxdata=new unsigned[size];
  memcpy(auxdata,data,sizeof(unsigned)*size);
  Sort(true,size,auxdata,nbits);
  delete[] auxdata;
}

//==============================================================================
/// Crea indice para ordenacion pero sin modificar los datos pasados.
//==============================================================================
void JRadixSort::MakeIndex(unsigned size,const unsigned long long *data,unsigned nbits){
  unsigned long long* auxdata=new unsigned long long[size];
  memcpy(auxdata,data,sizeof(unsigned long long)*size);
  Sort(true,size,auxdata,nbits);
  delete[] auxdata;
}

//==============================================================================
/// Ordena vector de datos en funcion del Index[] calculado previamente.
//==============================================================================
template<class T> void JRadixSort::TSortData(unsigned size,const T *data,T *result){
  const char met[]="TSortData";
  const int threads=omp_get_max_threads();
  if(!Index)RunException(met,"There is no index to sort data.");
  if(size!=Size)RunException(met,"The size of data is invalid.");
  T *res=result;
  if(data==res){//-Reserva vector auxiliar para la ordenacion.
    try{
      res=new T[size];
    }
    catch(const std::bad_alloc){
      RunException(met,"Cannot allocate the requested memory.");
    }
  }

  //-Reordena data[] en res[]
  if(!UseOmp || threads<2){//-Secuencial.
    for(unsigned c2=0;c2<Size;c2++)res[c2]=data[Index[c2]]; 
  }
  else{//-con OpenMP.
    const int nk=int(Size/OMPSIZE)+1;
    if(nk<0)RunException(met,"Number of values is invalid.");
    const int rk=int(Size%OMPSIZE);
    //printf(">> nk:%d  rk:%d\n",nk,rk);
    #ifdef _WITHOMP
      #pragma omp parallel for schedule (static)
    #endif
    for(int c=0;c<nk;c++){
      const unsigned c2ini=OMPSIZE*c;
      const unsigned c2fin=c2ini+(c+1<nk? OMPSIZE: rk);
      //printf(">> c2ini:%d  c2fin:%d\n",c2ini,c2fin);
      for(unsigned c2=c2ini;c2<c2fin;c2++)res[c2]=data[Index[c2]]; 
      //  printf("SortData> c2:%2u Index[c2]:%2u data[Index[c2]]:%g\n",c2,Index[c2],res[c2]);
    }
  }
  //-Coloca resultado y libera memoria.
  if(res!=result){
    memcpy(result,res,sizeof(T)*size);
    delete[] res;
  }
}
//==============================================================================
void JRadixSort::SortData(unsigned size,const byte *data,byte *result){ TSortData<byte>(size,data,result); }
void JRadixSort::SortData(unsigned size,const unsigned *data,unsigned *result){ TSortData<unsigned>(size,data,result); }
void JRadixSort::SortData(unsigned size,const float *data,float *result){ TSortData<float>(size,data,result); }
void JRadixSort::SortData(unsigned size,const double *data,double *result){ TSortData<double>(size,data,result); }
void JRadixSort::SortData(unsigned size,const tfloat3 *data,tfloat3 *result){ TSortData<tfloat3>(size,data,result); }
void JRadixSort::SortData(unsigned size,const tfloat4 *data,tfloat4 *result){ TSortData<tfloat4>(size,data,result); }
void JRadixSort::SortData(unsigned size,const tdouble3 *data,tdouble3 *result){ TSortData<tdouble3>(size,data,result); }
void JRadixSort::SortData(unsigned size,const tdouble2 *data,tdouble2 *result){ TSortData<tdouble2>(size,data,result); }

//==============================================================================
/// Comprueba ordenacion de datos.
//==============================================================================
void JRadixSort::DgCheckResult32()const{
  unsigned p=1;
  for(;p<Size&&InitData32[p-1]<=InitData32[p];p++);
  if(p!=Size)RunException("DgCheckResult32","El orden no es correcto");
}
//==============================================================================
/// Comprueba ordenacion de datos.
//==============================================================================
void JRadixSort::DgCheckResult64()const{
  unsigned p=1;
  for(;p<Size&&InitData64[p-1]<=InitData64[p];p++);
  if(p!=Size)RunException("DgCheckResult64","El orden no es correcto");
}











/*

//==============================================================================
/// Calcula Nbits para los datos indicados.
//==============================================================================
unsigned JRadixSort::CalcNbits(unsigned size,const unsigned *data){
  const char met[]="CalcNbits";
  const unsigned sizek=1024;
  const int nk=int(size/sizek);
  if(nk<0)RunException(met,"Numero de elementos no valido.");
  const int rk=int(size%sizek);
  //printf("CalcNbits> Size:%u  nk:%d  rk:%d  sizek:%u\n",size,nk,rk,sizek);
  //-Calcula maximo de nk bloques con varios hilos.
  unsigned *vmax=new unsigned[Threads*OMPSTRIDE];
  memset(vmax,0,sizeof(unsigned)*Threads*OMPSTRIDE);
#ifdef _WITHOMP
  #pragma omp parallel for schedule (static)
#endif
  for(int c=0;c<nk;c++){
    int th=omp_get_thread_num();
    unsigned mx=vmax[OMPSTRIDE*th];
    const unsigned c2fin=sizek*c+sizek;
    for(unsigned c2=c2fin-sizek;c2<c2fin;c2++)mx=max(mx,data[c2]);
    vmax[OMPSTRIDE*th]=mx;
  }
  //-Calcula reduce maximo de todos los hilos.
  unsigned mx=0;
  for(int t=0;t<Threads;t++)mx=max(mx,vmax[OMPSTRIDE*t]);
  delete[] vmax; vmax=NULL;
  //-Calcula maximo del resto de datos.
  const unsigned c2fin=sizek*nk+rk;
  for(unsigned c2=c2fin-rk;c2<c2fin;c2++)mx=max(mx,data[c2]);
  //-Calcula nbits para el valor maximo.
  return(BitsSize(mx));
}
//==============================================================================
/// Calcula Nbits para los datos indicados.
//==============================================================================
unsigned JRadixSort::CalcNbits(unsigned size,const unsigned long long *data){
  const char met[]="CalcNbits";
  const unsigned sizek=1024;
  const int nk=int(size/sizek);
  if(nk<0)RunException(met,"Number of values is invalid.");
  const int rk=int(size%sizek);
  //printf("CalcNbits> Size:%u  nk:%d  rk:%d  sizek:%u\n",size,nk,rk,sizek);
  //-Calcula maximo de nk bloques con varios hilos.
  unsigned long long *vmax=new unsigned long long[Threads*OMPSTRIDE];
  memset(vmax,0,sizeof(unsigned long long)*Threads*OMPSTRIDE);
  ThreadsConfig();
#ifdef _WITHOMP
  #pragma omp parallel for schedule (static)
#endif
  for(int c=0;c<nk;c++){
    int th=omp_get_thread_num();
    unsigned long long mx=vmax[OMPSTRIDE*th];
    const unsigned c2fin=sizek*c+sizek;
    for(unsigned c2=c2fin-sizek;c2<c2fin;c2++)mx=max(mx,data[c2]);
    vmax[OMPSTRIDE*th]=mx;
  }
  ThreadsRestore();
  //-Calcula reduce maximo de todos los hilos.
  unsigned long long mx=0;
  for(int t=0;t<Threads;t++)mx=max(mx,vmax[OMPSTRIDE*t]);
  delete[] vmax; vmax=NULL;
  //-Calcula maximo del resto de datos.
  const unsigned c2fin=sizek*nk+rk;
  for(unsigned c2=c2fin-rk;c2<c2fin;c2++)mx=max(mx,data[c2]);
  //-Calcula nbits para el valor maximo.
  return(BitsSize(mx));
}
//==============================================================================
/// Devuelve el numero de bits necesarios para codificar el valor indicado.
//==============================================================================
unsigned JRadixSort::BitsSize(unsigned v){
  unsigned nbits=1;
  for(;v>>nbits;nbits++);
  //sprintf(Cad,"++++> BitSize(%u):%u",v,nbits); Log->Print(Cad);
  return(nbits);
}
//==============================================================================
/// Devuelve el numero de bits necesarios para codificar el valor indicado.
//==============================================================================
unsigned JRadixSort::BitsSize(unsigned long long v){
  unsigned nbits=1;
  for(;v>>nbits;nbits++);
  //sprintf(Cad,"++++> BitSize(%u):%u",v,nbits); Log->Print(Cad);
  return(nbits);
}

*/


