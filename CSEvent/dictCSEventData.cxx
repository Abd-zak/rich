// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dictCSEventData

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "CSEventData.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_CSEventData(void *p = 0);
   static void *newArray_CSEventData(Long_t size, void *p);
   static void delete_CSEventData(void *p);
   static void deleteArray_CSEventData(void *p);
   static void destruct_CSEventData(void *p);
   static void streamer_CSEventData(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CSEventData*)
   {
      ::CSEventData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::CSEventData >(0);
      static ::ROOT::TGenericClassInfo 
         instance("CSEventData", ::CSEventData::Class_Version(), "CSEventData.h", 13,
                  typeid(::CSEventData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::CSEventData::Dictionary, isa_proxy, 16,
                  sizeof(::CSEventData) );
      instance.SetNew(&new_CSEventData);
      instance.SetNewArray(&newArray_CSEventData);
      instance.SetDelete(&delete_CSEventData);
      instance.SetDeleteArray(&deleteArray_CSEventData);
      instance.SetDestructor(&destruct_CSEventData);
      instance.SetStreamerFunc(&streamer_CSEventData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CSEventData*)
   {
      return GenerateInitInstanceLocal((::CSEventData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::CSEventData*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_CSHadronData(void *p = 0);
   static void *newArray_CSHadronData(Long_t size, void *p);
   static void delete_CSHadronData(void *p);
   static void deleteArray_CSHadronData(void *p);
   static void destruct_CSHadronData(void *p);
   static void streamer_CSHadronData(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CSHadronData*)
   {
      ::CSHadronData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::CSHadronData >(0);
      static ::ROOT::TGenericClassInfo 
         instance("CSHadronData", ::CSHadronData::Class_Version(), "CSEventData.h", 83,
                  typeid(::CSHadronData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::CSHadronData::Dictionary, isa_proxy, 16,
                  sizeof(::CSHadronData) );
      instance.SetNew(&new_CSHadronData);
      instance.SetNewArray(&newArray_CSHadronData);
      instance.SetDelete(&delete_CSHadronData);
      instance.SetDeleteArray(&deleteArray_CSHadronData);
      instance.SetDestructor(&destruct_CSHadronData);
      instance.SetStreamerFunc(&streamer_CSHadronData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CSHadronData*)
   {
      return GenerateInitInstanceLocal((::CSHadronData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::CSHadronData*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_CSResonanceData(void *p = 0);
   static void *newArray_CSResonanceData(Long_t size, void *p);
   static void delete_CSResonanceData(void *p);
   static void deleteArray_CSResonanceData(void *p);
   static void destruct_CSResonanceData(void *p);
   static void streamer_CSResonanceData(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CSResonanceData*)
   {
      ::CSResonanceData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::CSResonanceData >(0);
      static ::ROOT::TGenericClassInfo 
         instance("CSResonanceData", ::CSResonanceData::Class_Version(), "CSEventData.h", 141,
                  typeid(::CSResonanceData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::CSResonanceData::Dictionary, isa_proxy, 16,
                  sizeof(::CSResonanceData) );
      instance.SetNew(&new_CSResonanceData);
      instance.SetNewArray(&newArray_CSResonanceData);
      instance.SetDelete(&delete_CSResonanceData);
      instance.SetDeleteArray(&deleteArray_CSResonanceData);
      instance.SetDestructor(&destruct_CSResonanceData);
      instance.SetStreamerFunc(&streamer_CSResonanceData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CSResonanceData*)
   {
      return GenerateInitInstanceLocal((::CSResonanceData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::CSResonanceData*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr CSEventData::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *CSEventData::Class_Name()
{
   return "CSEventData";
}

//______________________________________________________________________________
const char *CSEventData::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CSEventData*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int CSEventData::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CSEventData*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *CSEventData::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CSEventData*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *CSEventData::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CSEventData*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr CSHadronData::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *CSHadronData::Class_Name()
{
   return "CSHadronData";
}

//______________________________________________________________________________
const char *CSHadronData::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CSHadronData*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int CSHadronData::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CSHadronData*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *CSHadronData::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CSHadronData*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *CSHadronData::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CSHadronData*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr CSResonanceData::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *CSResonanceData::Class_Name()
{
   return "CSResonanceData";
}

//______________________________________________________________________________
const char *CSResonanceData::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CSResonanceData*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int CSResonanceData::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CSResonanceData*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *CSResonanceData::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CSResonanceData*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *CSResonanceData::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CSResonanceData*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void CSEventData::Streamer(TBuffer &R__b)
{
   // Stream an object of class CSEventData.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      R__b >> runNo;
      R__b >> spillNo;
      R__b >> evtNo;
      R__b >> trigMask;
      R__b >> primaryPat;
      R__b >> spillOKPat;
      R__b >> Xp;
      R__b >> Yp;
      R__b >> Zp;
      R__b >> E0;
      R__b >> Q2;
      R__b >> y;
      R__b >> index;
      R__b >> prodIndex;
      R__b >> piThr;
      R__b >> dEK;
      R__b >> nOuts;
      R__b >> nKs;
      R__b >> nTrksRIt;
      R__b >> nTrksRIb;
      R__b >> KThr;
      R__b >> pThr;
      R__b.CheckByteCount(R__s, R__c, CSEventData::IsA());
   } else {
      R__c = R__b.WriteVersion(CSEventData::IsA(), kTRUE);
      R__b << runNo;
      R__b << spillNo;
      R__b << evtNo;
      R__b << trigMask;
      R__b << primaryPat;
      R__b << spillOKPat;
      R__b << Xp;
      R__b << Yp;
      R__b << Zp;
      R__b << E0;
      R__b << Q2;
      R__b << y;
      R__b << index;
      R__b << prodIndex;
      R__b << piThr;
      R__b << dEK;
      R__b << nOuts;
      R__b << nKs;
      R__b << nTrksRIt;
      R__b << nTrksRIb;
      R__b << KThr;
      R__b << pThr;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_CSEventData(void *p) {
      return  p ? new(p) ::CSEventData : new ::CSEventData;
   }
   static void *newArray_CSEventData(Long_t nElements, void *p) {
      return p ? new(p) ::CSEventData[nElements] : new ::CSEventData[nElements];
   }
   // Wrapper around operator delete
   static void delete_CSEventData(void *p) {
      delete ((::CSEventData*)p);
   }
   static void deleteArray_CSEventData(void *p) {
      delete [] ((::CSEventData*)p);
   }
   static void destruct_CSEventData(void *p) {
      typedef ::CSEventData current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_CSEventData(TBuffer &buf, void *obj) {
      ((::CSEventData*)obj)->::CSEventData::Streamer(buf);
   }
} // end of namespace ROOT for class ::CSEventData

//______________________________________________________________________________
void CSHadronData::Streamer(TBuffer &R__b)
{
   // Stream an object of class CSHadronData.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      R__b >> Px;
      R__b >> Py;
      R__b >> Pz;
      R__b >> qP;
      R__b >> XX0;
      R__b >> ZFirst;
      R__b >> ZLast;
      R__b >> chi2;
      R__b >> ECAL;
      R__b >> phiR;
      R__b.ReadStaticArray((float*)LH);
      R__b.ReadStaticArray((float*)dLHdI);
      R__b >> thC;
      R__b >> MCpid;
      R__b >> hasR;
      R__b >> XR;
      R__b >> YR;
      R__b >> tgXR;
      R__b >> tgYR;
      R__b >> P;
      R__b >> q;
      R__b >> thR;
      R__b.CheckByteCount(R__s, R__c, CSHadronData::IsA());
   } else {
      R__c = R__b.WriteVersion(CSHadronData::IsA(), kTRUE);
      R__b << Px;
      R__b << Py;
      R__b << Pz;
      R__b << qP;
      R__b << XX0;
      R__b << ZFirst;
      R__b << ZLast;
      R__b << chi2;
      R__b << ECAL;
      R__b << phiR;
      R__b.WriteArray(LH, 6);
      R__b.WriteArray(dLHdI, 3);
      R__b << thC;
      R__b << MCpid;
      R__b << hasR;
      R__b << XR;
      R__b << YR;
      R__b << tgXR;
      R__b << tgYR;
      R__b << P;
      R__b << q;
      R__b << thR;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_CSHadronData(void *p) {
      return  p ? new(p) ::CSHadronData : new ::CSHadronData;
   }
   static void *newArray_CSHadronData(Long_t nElements, void *p) {
      return p ? new(p) ::CSHadronData[nElements] : new ::CSHadronData[nElements];
   }
   // Wrapper around operator delete
   static void delete_CSHadronData(void *p) {
      delete ((::CSHadronData*)p);
   }
   static void deleteArray_CSHadronData(void *p) {
      delete [] ((::CSHadronData*)p);
   }
   static void destruct_CSHadronData(void *p) {
      typedef ::CSHadronData current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_CSHadronData(TBuffer &buf, void *obj) {
      ((::CSHadronData*)obj)->::CSHadronData::Streamer(buf);
   }
} // end of namespace ROOT for class ::CSHadronData

//______________________________________________________________________________
void CSResonanceData::Streamer(TBuffer &R__b)
{
   // Stream an object of class CSResonanceData.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      R__b >> phiPat;
      R__b >> K0Pat;
      R__b >> LambdaPat;
      R__b >> h1;
      R__b >> h2;
      R__b >> m;
      R__b >> Xs;
      R__b >> Ys;
      R__b >> Zs;
      R__b >> D;
      R__b >> dD;
      R__b >> cth;
      R__b >> pT;
      R__b >> alpha;
      R__b >> chi2;
      R__b.CheckByteCount(R__s, R__c, CSResonanceData::IsA());
   } else {
      R__c = R__b.WriteVersion(CSResonanceData::IsA(), kTRUE);
      R__b << phiPat;
      R__b << K0Pat;
      R__b << LambdaPat;
      R__b << h1;
      R__b << h2;
      R__b << m;
      R__b << Xs;
      R__b << Ys;
      R__b << Zs;
      R__b << D;
      R__b << dD;
      R__b << cth;
      R__b << pT;
      R__b << alpha;
      R__b << chi2;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_CSResonanceData(void *p) {
      return  p ? new(p) ::CSResonanceData : new ::CSResonanceData;
   }
   static void *newArray_CSResonanceData(Long_t nElements, void *p) {
      return p ? new(p) ::CSResonanceData[nElements] : new ::CSResonanceData[nElements];
   }
   // Wrapper around operator delete
   static void delete_CSResonanceData(void *p) {
      delete ((::CSResonanceData*)p);
   }
   static void deleteArray_CSResonanceData(void *p) {
      delete [] ((::CSResonanceData*)p);
   }
   static void destruct_CSResonanceData(void *p) {
      typedef ::CSResonanceData current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_CSResonanceData(TBuffer &buf, void *obj) {
      ((::CSResonanceData*)obj)->::CSResonanceData::Streamer(buf);
   }
} // end of namespace ROOT for class ::CSResonanceData

namespace ROOT {
   static TClass *vectorlECSResonanceDatagR_Dictionary();
   static void vectorlECSResonanceDatagR_TClassManip(TClass*);
   static void *new_vectorlECSResonanceDatagR(void *p = 0);
   static void *newArray_vectorlECSResonanceDatagR(Long_t size, void *p);
   static void delete_vectorlECSResonanceDatagR(void *p);
   static void deleteArray_vectorlECSResonanceDatagR(void *p);
   static void destruct_vectorlECSResonanceDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<CSResonanceData>*)
   {
      vector<CSResonanceData> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<CSResonanceData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<CSResonanceData>", -2, "vector", 214,
                  typeid(vector<CSResonanceData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlECSResonanceDatagR_Dictionary, isa_proxy, 4,
                  sizeof(vector<CSResonanceData>) );
      instance.SetNew(&new_vectorlECSResonanceDatagR);
      instance.SetNewArray(&newArray_vectorlECSResonanceDatagR);
      instance.SetDelete(&delete_vectorlECSResonanceDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlECSResonanceDatagR);
      instance.SetDestructor(&destruct_vectorlECSResonanceDatagR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<CSResonanceData> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<CSResonanceData>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlECSResonanceDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<CSResonanceData>*)0x0)->GetClass();
      vectorlECSResonanceDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlECSResonanceDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlECSResonanceDatagR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<CSResonanceData> : new vector<CSResonanceData>;
   }
   static void *newArray_vectorlECSResonanceDatagR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<CSResonanceData>[nElements] : new vector<CSResonanceData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlECSResonanceDatagR(void *p) {
      delete ((vector<CSResonanceData>*)p);
   }
   static void deleteArray_vectorlECSResonanceDatagR(void *p) {
      delete [] ((vector<CSResonanceData>*)p);
   }
   static void destruct_vectorlECSResonanceDatagR(void *p) {
      typedef vector<CSResonanceData> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<CSResonanceData>

namespace ROOT {
   static TClass *vectorlECSHadronDatagR_Dictionary();
   static void vectorlECSHadronDatagR_TClassManip(TClass*);
   static void *new_vectorlECSHadronDatagR(void *p = 0);
   static void *newArray_vectorlECSHadronDatagR(Long_t size, void *p);
   static void delete_vectorlECSHadronDatagR(void *p);
   static void deleteArray_vectorlECSHadronDatagR(void *p);
   static void destruct_vectorlECSHadronDatagR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<CSHadronData>*)
   {
      vector<CSHadronData> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<CSHadronData>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<CSHadronData>", -2, "vector", 214,
                  typeid(vector<CSHadronData>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlECSHadronDatagR_Dictionary, isa_proxy, 4,
                  sizeof(vector<CSHadronData>) );
      instance.SetNew(&new_vectorlECSHadronDatagR);
      instance.SetNewArray(&newArray_vectorlECSHadronDatagR);
      instance.SetDelete(&delete_vectorlECSHadronDatagR);
      instance.SetDeleteArray(&deleteArray_vectorlECSHadronDatagR);
      instance.SetDestructor(&destruct_vectorlECSHadronDatagR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<CSHadronData> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<CSHadronData>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlECSHadronDatagR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<CSHadronData>*)0x0)->GetClass();
      vectorlECSHadronDatagR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlECSHadronDatagR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlECSHadronDatagR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<CSHadronData> : new vector<CSHadronData>;
   }
   static void *newArray_vectorlECSHadronDatagR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<CSHadronData>[nElements] : new vector<CSHadronData>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlECSHadronDatagR(void *p) {
      delete ((vector<CSHadronData>*)p);
   }
   static void deleteArray_vectorlECSHadronDatagR(void *p) {
      delete [] ((vector<CSHadronData>*)p);
   }
   static void destruct_vectorlECSHadronDatagR(void *p) {
      typedef vector<CSHadronData> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<CSHadronData>

namespace {
  void TriggerDictionaryInitialization_dictCSEventData_Impl() {
    static const char* headers[] = {
"CSEventData.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.14.04-4d676/x86_64-centos7-gcc62-opt/include",
"/pbs/home/a/azakaria/rich_matrices/CSEvent/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "dictCSEventData dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
struct __attribute__((annotate(R"ATTRDUMP(Must be the last item before the closing '};')ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Must be the last item before the closing '};')ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Must be the last item before the closing '};')ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Must be the last item before the closing '};')ATTRDUMP"))) __attribute__((annotate("$clingAutoload$CSEventData.h")))  CSResonanceData;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
struct __attribute__((annotate(R"ATTRDUMP(Must be the last item before the closing '};')ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Must be the last item before the closing '};')ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Must be the last item before the closing '};')ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Must be the last item before the closing '};')ATTRDUMP"))) __attribute__((annotate("$clingAutoload$CSEventData.h")))  CSHadronData;
struct __attribute__((annotate(R"ATTRDUMP(Must be the last item before the closing '};')ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Must be the last item before the closing '};')ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Must be the last item before the closing '};')ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Must be the last item before the closing '};')ATTRDUMP"))) __attribute__((annotate("$clingAutoload$CSEventData.h")))  CSEventData;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "dictCSEventData dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "CSEventData.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"CSEventData", payloadCode, "@",
"CSHadronData", payloadCode, "@",
"CSResonanceData", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("dictCSEventData",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_dictCSEventData_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_dictCSEventData_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_dictCSEventData() {
  TriggerDictionaryInitialization_dictCSEventData_Impl();
}
