// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME RooRarFitCint

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
#include "rarStrParser.hh"
#include "RooLass.hh"
#include "rarTwoGauss.hh"
#include "rarMultPdf.hh"
#include "RooOsipDisc.hh"
#include "rarStep.hh"
#include "rarProd.hh"
#include "rarNLL.hh"
#include "rarUniform.hh"
#include "rarConfig.hh"
#include "rarOsipDisc.hh"
#include "RooBallack.hh"
#include "rarLass.hh"
#include "rarAdd.hh"
#include "rarTriGauss.hh"
#include "rarNovosibirsk.hh"
#include "rarHistPdf.hh"
#include "rarSPlot.hh"
#include "RooThreshold.hh"
#include "rarGaussian.hh"
#include "RooCruijff.hh"
#include "rarBasePdf.hh"
#include "rarDatasets.hh"
#include "rarBinned.hh"
#include "rarBallack.hh"
#include "rarRelBreitWigner.hh"
#include "rarBifurGauss.hh"
#include "rarUsrPdf.hh"
#include "rarDatasetDef.hh"
#include "RooGounarisSakurai.hh"
#include "rarCompBase.hh"
#include "rarCBShape.hh"
#include "rarGeneric.hh"
#include "RooBinnedPdf.hh"
#include "rarCruijff.hh"
#include "rarMinuit.hh"
#include "rarSimPdf.hh"
#include "rarFlatte.hh"
#include "rarVoigtian.hh"
#include "rarThreshold.hh"
#include "rarGaussModel.hh"
#include "rarExp.hh"
#include "rarToyList.hh"
#include "rarArgusBG.hh"
#include "RooRelBreitWigner.hh"
#include "rarGounarisSakurai.hh"
#include "rarMLFitter.hh"
#include "rarKeys.hh"
#include "rarMLPdf.hh"
#include "RooFlatte.hh"
#include "rarDecay.hh"
#include "rarPoly.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_rarStrParser(void *p = 0);
   static void *newArray_rarStrParser(Long_t size, void *p);
   static void delete_rarStrParser(void *p);
   static void deleteArray_rarStrParser(void *p);
   static void destruct_rarStrParser(void *p);
   static void streamer_rarStrParser(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarStrParser*)
   {
      ::rarStrParser *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarStrParser >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarStrParser", ::rarStrParser::Class_Version(), "rarStrParser.hh", 27,
                  typeid(::rarStrParser), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarStrParser::Dictionary, isa_proxy, 16,
                  sizeof(::rarStrParser) );
      instance.SetNew(&new_rarStrParser);
      instance.SetNewArray(&newArray_rarStrParser);
      instance.SetDelete(&delete_rarStrParser);
      instance.SetDeleteArray(&deleteArray_rarStrParser);
      instance.SetDestructor(&destruct_rarStrParser);
      instance.SetStreamerFunc(&streamer_rarStrParser);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarStrParser*)
   {
      return GenerateInitInstanceLocal((::rarStrParser*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarStrParser*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_RooLass(void *p);
   static void deleteArray_RooLass(void *p);
   static void destruct_RooLass(void *p);
   static void streamer_RooLass(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooLass*)
   {
      ::RooLass *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooLass >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooLass", ::RooLass::Class_Version(), "RooLass.hh", 21,
                  typeid(::RooLass), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooLass::Dictionary, isa_proxy, 16,
                  sizeof(::RooLass) );
      instance.SetDelete(&delete_RooLass);
      instance.SetDeleteArray(&deleteArray_RooLass);
      instance.SetDestructor(&destruct_RooLass);
      instance.SetStreamerFunc(&streamer_RooLass);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooLass*)
   {
      return GenerateInitInstanceLocal((::RooLass*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RooLass*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarConfig(void *p = 0);
   static void *newArray_rarConfig(Long_t size, void *p);
   static void delete_rarConfig(void *p);
   static void deleteArray_rarConfig(void *p);
   static void destruct_rarConfig(void *p);
   static void streamer_rarConfig(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarConfig*)
   {
      ::rarConfig *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarConfig >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarConfig", ::rarConfig::Class_Version(), "rarConfig.hh", 27,
                  typeid(::rarConfig), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarConfig::Dictionary, isa_proxy, 16,
                  sizeof(::rarConfig) );
      instance.SetNew(&new_rarConfig);
      instance.SetNewArray(&newArray_rarConfig);
      instance.SetDelete(&delete_rarConfig);
      instance.SetDeleteArray(&deleteArray_rarConfig);
      instance.SetDestructor(&destruct_rarConfig);
      instance.SetStreamerFunc(&streamer_rarConfig);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarConfig*)
   {
      return GenerateInitInstanceLocal((::rarConfig*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarConfig*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarDatasetDef(void *p = 0);
   static void *newArray_rarDatasetDef(Long_t size, void *p);
   static void delete_rarDatasetDef(void *p);
   static void deleteArray_rarDatasetDef(void *p);
   static void destruct_rarDatasetDef(void *p);
   static void streamer_rarDatasetDef(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarDatasetDef*)
   {
      ::rarDatasetDef *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarDatasetDef >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarDatasetDef", ::rarDatasetDef::Class_Version(), "rarDatasetDef.hh", 47,
                  typeid(::rarDatasetDef), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarDatasetDef::Dictionary, isa_proxy, 16,
                  sizeof(::rarDatasetDef) );
      instance.SetNew(&new_rarDatasetDef);
      instance.SetNewArray(&newArray_rarDatasetDef);
      instance.SetDelete(&delete_rarDatasetDef);
      instance.SetDeleteArray(&deleteArray_rarDatasetDef);
      instance.SetDestructor(&destruct_rarDatasetDef);
      instance.SetStreamerFunc(&streamer_rarDatasetDef);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarDatasetDef*)
   {
      return GenerateInitInstanceLocal((::rarDatasetDef*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarDatasetDef*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarDatasets(void *p = 0);
   static void *newArray_rarDatasets(Long_t size, void *p);
   static void delete_rarDatasets(void *p);
   static void deleteArray_rarDatasets(void *p);
   static void destruct_rarDatasets(void *p);
   static void streamer_rarDatasets(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarDatasets*)
   {
      ::rarDatasets *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarDatasets >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarDatasets", ::rarDatasets::Class_Version(), "rarDatasets.hh", 27,
                  typeid(::rarDatasets), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarDatasets::Dictionary, isa_proxy, 16,
                  sizeof(::rarDatasets) );
      instance.SetNew(&new_rarDatasets);
      instance.SetNewArray(&newArray_rarDatasets);
      instance.SetDelete(&delete_rarDatasets);
      instance.SetDeleteArray(&deleteArray_rarDatasets);
      instance.SetDestructor(&destruct_rarDatasets);
      instance.SetStreamerFunc(&streamer_rarDatasets);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarDatasets*)
   {
      return GenerateInitInstanceLocal((::rarDatasets*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarDatasets*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarBasePdf(void *p = 0);
   static void *newArray_rarBasePdf(Long_t size, void *p);
   static void delete_rarBasePdf(void *p);
   static void deleteArray_rarBasePdf(void *p);
   static void destruct_rarBasePdf(void *p);
   static void streamer_rarBasePdf(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarBasePdf*)
   {
      ::rarBasePdf *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarBasePdf >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarBasePdf", ::rarBasePdf::Class_Version(), "rarBasePdf.hh", 39,
                  typeid(::rarBasePdf), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarBasePdf::Dictionary, isa_proxy, 16,
                  sizeof(::rarBasePdf) );
      instance.SetNew(&new_rarBasePdf);
      instance.SetNewArray(&newArray_rarBasePdf);
      instance.SetDelete(&delete_rarBasePdf);
      instance.SetDeleteArray(&deleteArray_rarBasePdf);
      instance.SetDestructor(&destruct_rarBasePdf);
      instance.SetStreamerFunc(&streamer_rarBasePdf);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarBasePdf*)
   {
      return GenerateInitInstanceLocal((::rarBasePdf*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarBasePdf*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarTwoGauss(void *p = 0);
   static void *newArray_rarTwoGauss(Long_t size, void *p);
   static void delete_rarTwoGauss(void *p);
   static void deleteArray_rarTwoGauss(void *p);
   static void destruct_rarTwoGauss(void *p);
   static void streamer_rarTwoGauss(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarTwoGauss*)
   {
      ::rarTwoGauss *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarTwoGauss >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarTwoGauss", ::rarTwoGauss::Class_Version(), "rarTwoGauss.hh", 24,
                  typeid(::rarTwoGauss), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarTwoGauss::Dictionary, isa_proxy, 16,
                  sizeof(::rarTwoGauss) );
      instance.SetNew(&new_rarTwoGauss);
      instance.SetNewArray(&newArray_rarTwoGauss);
      instance.SetDelete(&delete_rarTwoGauss);
      instance.SetDeleteArray(&deleteArray_rarTwoGauss);
      instance.SetDestructor(&destruct_rarTwoGauss);
      instance.SetStreamerFunc(&streamer_rarTwoGauss);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarTwoGauss*)
   {
      return GenerateInitInstanceLocal((::rarTwoGauss*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarTwoGauss*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarCompBase(void *p = 0);
   static void *newArray_rarCompBase(Long_t size, void *p);
   static void delete_rarCompBase(void *p);
   static void deleteArray_rarCompBase(void *p);
   static void destruct_rarCompBase(void *p);
   static void streamer_rarCompBase(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarCompBase*)
   {
      ::rarCompBase *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarCompBase >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarCompBase", ::rarCompBase::Class_Version(), "rarCompBase.hh", 35,
                  typeid(::rarCompBase), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarCompBase::Dictionary, isa_proxy, 16,
                  sizeof(::rarCompBase) );
      instance.SetNew(&new_rarCompBase);
      instance.SetNewArray(&newArray_rarCompBase);
      instance.SetDelete(&delete_rarCompBase);
      instance.SetDeleteArray(&deleteArray_rarCompBase);
      instance.SetDestructor(&destruct_rarCompBase);
      instance.SetStreamerFunc(&streamer_rarCompBase);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarCompBase*)
   {
      return GenerateInitInstanceLocal((::rarCompBase*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarCompBase*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarMultPdf(void *p = 0);
   static void *newArray_rarMultPdf(Long_t size, void *p);
   static void delete_rarMultPdf(void *p);
   static void deleteArray_rarMultPdf(void *p);
   static void destruct_rarMultPdf(void *p);
   static void streamer_rarMultPdf(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarMultPdf*)
   {
      ::rarMultPdf *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarMultPdf >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarMultPdf", ::rarMultPdf::Class_Version(), "rarMultPdf.hh", 26,
                  typeid(::rarMultPdf), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarMultPdf::Dictionary, isa_proxy, 16,
                  sizeof(::rarMultPdf) );
      instance.SetNew(&new_rarMultPdf);
      instance.SetNewArray(&newArray_rarMultPdf);
      instance.SetDelete(&delete_rarMultPdf);
      instance.SetDeleteArray(&deleteArray_rarMultPdf);
      instance.SetDestructor(&destruct_rarMultPdf);
      instance.SetStreamerFunc(&streamer_rarMultPdf);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarMultPdf*)
   {
      return GenerateInitInstanceLocal((::rarMultPdf*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarMultPdf*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_RooOsipDisc(void *p);
   static void deleteArray_RooOsipDisc(void *p);
   static void destruct_RooOsipDisc(void *p);
   static void streamer_RooOsipDisc(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooOsipDisc*)
   {
      ::RooOsipDisc *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooOsipDisc >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooOsipDisc", ::RooOsipDisc::Class_Version(), "RooOsipDisc.hh", 24,
                  typeid(::RooOsipDisc), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooOsipDisc::Dictionary, isa_proxy, 16,
                  sizeof(::RooOsipDisc) );
      instance.SetDelete(&delete_RooOsipDisc);
      instance.SetDeleteArray(&deleteArray_RooOsipDisc);
      instance.SetDestructor(&destruct_RooOsipDisc);
      instance.SetStreamerFunc(&streamer_RooOsipDisc);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooOsipDisc*)
   {
      return GenerateInitInstanceLocal((::RooOsipDisc*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RooOsipDisc*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarStep(void *p = 0);
   static void *newArray_rarStep(Long_t size, void *p);
   static void delete_rarStep(void *p);
   static void deleteArray_rarStep(void *p);
   static void destruct_rarStep(void *p);
   static void streamer_rarStep(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarStep*)
   {
      ::rarStep *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarStep >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarStep", ::rarStep::Class_Version(), "rarStep.hh", 42,
                  typeid(::rarStep), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarStep::Dictionary, isa_proxy, 16,
                  sizeof(::rarStep) );
      instance.SetNew(&new_rarStep);
      instance.SetNewArray(&newArray_rarStep);
      instance.SetDelete(&delete_rarStep);
      instance.SetDeleteArray(&deleteArray_rarStep);
      instance.SetDestructor(&destruct_rarStep);
      instance.SetStreamerFunc(&streamer_rarStep);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarStep*)
   {
      return GenerateInitInstanceLocal((::rarStep*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarStep*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarProd(void *p = 0);
   static void *newArray_rarProd(Long_t size, void *p);
   static void delete_rarProd(void *p);
   static void deleteArray_rarProd(void *p);
   static void destruct_rarProd(void *p);
   static void streamer_rarProd(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarProd*)
   {
      ::rarProd *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarProd >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarProd", ::rarProd::Class_Version(), "rarProd.hh", 26,
                  typeid(::rarProd), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarProd::Dictionary, isa_proxy, 16,
                  sizeof(::rarProd) );
      instance.SetNew(&new_rarProd);
      instance.SetNewArray(&newArray_rarProd);
      instance.SetDelete(&delete_rarProd);
      instance.SetDeleteArray(&deleteArray_rarProd);
      instance.SetDestructor(&destruct_rarProd);
      instance.SetStreamerFunc(&streamer_rarProd);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarProd*)
   {
      return GenerateInitInstanceLocal((::rarProd*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarProd*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarNLL(void *p = 0);
   static void *newArray_rarNLL(Long_t size, void *p);
   static void delete_rarNLL(void *p);
   static void deleteArray_rarNLL(void *p);
   static void destruct_rarNLL(void *p);
   static void streamer_rarNLL(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarNLL*)
   {
      ::rarNLL *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarNLL >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarNLL", ::rarNLL::Class_Version(), "rarNLL.hh", 25,
                  typeid(::rarNLL), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarNLL::Dictionary, isa_proxy, 16,
                  sizeof(::rarNLL) );
      instance.SetNew(&new_rarNLL);
      instance.SetNewArray(&newArray_rarNLL);
      instance.SetDelete(&delete_rarNLL);
      instance.SetDeleteArray(&deleteArray_rarNLL);
      instance.SetDestructor(&destruct_rarNLL);
      instance.SetStreamerFunc(&streamer_rarNLL);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarNLL*)
   {
      return GenerateInitInstanceLocal((::rarNLL*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarNLL*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarUniform(void *p = 0);
   static void *newArray_rarUniform(Long_t size, void *p);
   static void delete_rarUniform(void *p);
   static void deleteArray_rarUniform(void *p);
   static void destruct_rarUniform(void *p);
   static void streamer_rarUniform(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarUniform*)
   {
      ::rarUniform *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarUniform >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarUniform", ::rarUniform::Class_Version(), "rarUniform.hh", 31,
                  typeid(::rarUniform), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarUniform::Dictionary, isa_proxy, 16,
                  sizeof(::rarUniform) );
      instance.SetNew(&new_rarUniform);
      instance.SetNewArray(&newArray_rarUniform);
      instance.SetDelete(&delete_rarUniform);
      instance.SetDeleteArray(&deleteArray_rarUniform);
      instance.SetDestructor(&destruct_rarUniform);
      instance.SetStreamerFunc(&streamer_rarUniform);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarUniform*)
   {
      return GenerateInitInstanceLocal((::rarUniform*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarUniform*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarOsipDisc(void *p = 0);
   static void *newArray_rarOsipDisc(Long_t size, void *p);
   static void delete_rarOsipDisc(void *p);
   static void deleteArray_rarOsipDisc(void *p);
   static void destruct_rarOsipDisc(void *p);
   static void streamer_rarOsipDisc(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarOsipDisc*)
   {
      ::rarOsipDisc *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarOsipDisc >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarOsipDisc", ::rarOsipDisc::Class_Version(), "rarOsipDisc.hh", 38,
                  typeid(::rarOsipDisc), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarOsipDisc::Dictionary, isa_proxy, 16,
                  sizeof(::rarOsipDisc) );
      instance.SetNew(&new_rarOsipDisc);
      instance.SetNewArray(&newArray_rarOsipDisc);
      instance.SetDelete(&delete_rarOsipDisc);
      instance.SetDeleteArray(&deleteArray_rarOsipDisc);
      instance.SetDestructor(&destruct_rarOsipDisc);
      instance.SetStreamerFunc(&streamer_rarOsipDisc);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarOsipDisc*)
   {
      return GenerateInitInstanceLocal((::rarOsipDisc*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarOsipDisc*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_RooBallack(void *p);
   static void deleteArray_RooBallack(void *p);
   static void destruct_RooBallack(void *p);
   static void streamer_RooBallack(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooBallack*)
   {
      ::RooBallack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooBallack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooBallack", ::RooBallack::Class_Version(), "RooBallack.hh", 19,
                  typeid(::RooBallack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooBallack::Dictionary, isa_proxy, 16,
                  sizeof(::RooBallack) );
      instance.SetDelete(&delete_RooBallack);
      instance.SetDeleteArray(&deleteArray_RooBallack);
      instance.SetDestructor(&destruct_RooBallack);
      instance.SetStreamerFunc(&streamer_RooBallack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooBallack*)
   {
      return GenerateInitInstanceLocal((::RooBallack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RooBallack*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarLass(void *p = 0);
   static void *newArray_rarLass(Long_t size, void *p);
   static void delete_rarLass(void *p);
   static void deleteArray_rarLass(void *p);
   static void destruct_rarLass(void *p);
   static void streamer_rarLass(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarLass*)
   {
      ::rarLass *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarLass >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarLass", ::rarLass::Class_Version(), "rarLass.hh", 24,
                  typeid(::rarLass), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarLass::Dictionary, isa_proxy, 16,
                  sizeof(::rarLass) );
      instance.SetNew(&new_rarLass);
      instance.SetNewArray(&newArray_rarLass);
      instance.SetDelete(&delete_rarLass);
      instance.SetDeleteArray(&deleteArray_rarLass);
      instance.SetDestructor(&destruct_rarLass);
      instance.SetStreamerFunc(&streamer_rarLass);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarLass*)
   {
      return GenerateInitInstanceLocal((::rarLass*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarLass*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarAdd(void *p = 0);
   static void *newArray_rarAdd(Long_t size, void *p);
   static void delete_rarAdd(void *p);
   static void deleteArray_rarAdd(void *p);
   static void destruct_rarAdd(void *p);
   static void streamer_rarAdd(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarAdd*)
   {
      ::rarAdd *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarAdd >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarAdd", ::rarAdd::Class_Version(), "rarAdd.hh", 28,
                  typeid(::rarAdd), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarAdd::Dictionary, isa_proxy, 16,
                  sizeof(::rarAdd) );
      instance.SetNew(&new_rarAdd);
      instance.SetNewArray(&newArray_rarAdd);
      instance.SetDelete(&delete_rarAdd);
      instance.SetDeleteArray(&deleteArray_rarAdd);
      instance.SetDestructor(&destruct_rarAdd);
      instance.SetStreamerFunc(&streamer_rarAdd);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarAdd*)
   {
      return GenerateInitInstanceLocal((::rarAdd*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarAdd*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarTriGauss(void *p = 0);
   static void *newArray_rarTriGauss(Long_t size, void *p);
   static void delete_rarTriGauss(void *p);
   static void deleteArray_rarTriGauss(void *p);
   static void destruct_rarTriGauss(void *p);
   static void streamer_rarTriGauss(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarTriGauss*)
   {
      ::rarTriGauss *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarTriGauss >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarTriGauss", ::rarTriGauss::Class_Version(), "rarTriGauss.hh", 26,
                  typeid(::rarTriGauss), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarTriGauss::Dictionary, isa_proxy, 16,
                  sizeof(::rarTriGauss) );
      instance.SetNew(&new_rarTriGauss);
      instance.SetNewArray(&newArray_rarTriGauss);
      instance.SetDelete(&delete_rarTriGauss);
      instance.SetDeleteArray(&deleteArray_rarTriGauss);
      instance.SetDestructor(&destruct_rarTriGauss);
      instance.SetStreamerFunc(&streamer_rarTriGauss);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarTriGauss*)
   {
      return GenerateInitInstanceLocal((::rarTriGauss*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarTriGauss*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarNovosibirsk(void *p = 0);
   static void *newArray_rarNovosibirsk(Long_t size, void *p);
   static void delete_rarNovosibirsk(void *p);
   static void deleteArray_rarNovosibirsk(void *p);
   static void destruct_rarNovosibirsk(void *p);
   static void streamer_rarNovosibirsk(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarNovosibirsk*)
   {
      ::rarNovosibirsk *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarNovosibirsk >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarNovosibirsk", ::rarNovosibirsk::Class_Version(), "rarNovosibirsk.hh", 36,
                  typeid(::rarNovosibirsk), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarNovosibirsk::Dictionary, isa_proxy, 16,
                  sizeof(::rarNovosibirsk) );
      instance.SetNew(&new_rarNovosibirsk);
      instance.SetNewArray(&newArray_rarNovosibirsk);
      instance.SetDelete(&delete_rarNovosibirsk);
      instance.SetDeleteArray(&deleteArray_rarNovosibirsk);
      instance.SetDestructor(&destruct_rarNovosibirsk);
      instance.SetStreamerFunc(&streamer_rarNovosibirsk);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarNovosibirsk*)
   {
      return GenerateInitInstanceLocal((::rarNovosibirsk*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarNovosibirsk*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarHistPdf(void *p = 0);
   static void *newArray_rarHistPdf(Long_t size, void *p);
   static void delete_rarHistPdf(void *p);
   static void deleteArray_rarHistPdf(void *p);
   static void destruct_rarHistPdf(void *p);
   static void streamer_rarHistPdf(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarHistPdf*)
   {
      ::rarHistPdf *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarHistPdf >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarHistPdf", ::rarHistPdf::Class_Version(), "rarHistPdf.hh", 28,
                  typeid(::rarHistPdf), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarHistPdf::Dictionary, isa_proxy, 16,
                  sizeof(::rarHistPdf) );
      instance.SetNew(&new_rarHistPdf);
      instance.SetNewArray(&newArray_rarHistPdf);
      instance.SetDelete(&delete_rarHistPdf);
      instance.SetDeleteArray(&deleteArray_rarHistPdf);
      instance.SetDestructor(&destruct_rarHistPdf);
      instance.SetStreamerFunc(&streamer_rarHistPdf);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarHistPdf*)
   {
      return GenerateInitInstanceLocal((::rarHistPdf*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarHistPdf*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarSPlot(void *p = 0);
   static void *newArray_rarSPlot(Long_t size, void *p);
   static void delete_rarSPlot(void *p);
   static void deleteArray_rarSPlot(void *p);
   static void destruct_rarSPlot(void *p);
   static void streamer_rarSPlot(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarSPlot*)
   {
      ::rarSPlot *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarSPlot >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarSPlot", ::rarSPlot::Class_Version(), "rarSPlot.hh", 45,
                  typeid(::rarSPlot), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarSPlot::Dictionary, isa_proxy, 16,
                  sizeof(::rarSPlot) );
      instance.SetNew(&new_rarSPlot);
      instance.SetNewArray(&newArray_rarSPlot);
      instance.SetDelete(&delete_rarSPlot);
      instance.SetDeleteArray(&deleteArray_rarSPlot);
      instance.SetDestructor(&destruct_rarSPlot);
      instance.SetStreamerFunc(&streamer_rarSPlot);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarSPlot*)
   {
      return GenerateInitInstanceLocal((::rarSPlot*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarSPlot*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_RooThreshold(void *p);
   static void deleteArray_RooThreshold(void *p);
   static void destruct_RooThreshold(void *p);
   static void streamer_RooThreshold(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooThreshold*)
   {
      ::RooThreshold *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooThreshold >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooThreshold", ::RooThreshold::Class_Version(), "RooThreshold.hh", 20,
                  typeid(::RooThreshold), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooThreshold::Dictionary, isa_proxy, 16,
                  sizeof(::RooThreshold) );
      instance.SetDelete(&delete_RooThreshold);
      instance.SetDeleteArray(&deleteArray_RooThreshold);
      instance.SetDestructor(&destruct_RooThreshold);
      instance.SetStreamerFunc(&streamer_RooThreshold);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooThreshold*)
   {
      return GenerateInitInstanceLocal((::RooThreshold*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RooThreshold*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarGaussian(void *p = 0);
   static void *newArray_rarGaussian(Long_t size, void *p);
   static void delete_rarGaussian(void *p);
   static void deleteArray_rarGaussian(void *p);
   static void destruct_rarGaussian(void *p);
   static void streamer_rarGaussian(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarGaussian*)
   {
      ::rarGaussian *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarGaussian >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarGaussian", ::rarGaussian::Class_Version(), "rarGaussian.hh", 38,
                  typeid(::rarGaussian), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarGaussian::Dictionary, isa_proxy, 16,
                  sizeof(::rarGaussian) );
      instance.SetNew(&new_rarGaussian);
      instance.SetNewArray(&newArray_rarGaussian);
      instance.SetDelete(&delete_rarGaussian);
      instance.SetDeleteArray(&deleteArray_rarGaussian);
      instance.SetDestructor(&destruct_rarGaussian);
      instance.SetStreamerFunc(&streamer_rarGaussian);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarGaussian*)
   {
      return GenerateInitInstanceLocal((::rarGaussian*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarGaussian*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_RooCruijff(void *p);
   static void deleteArray_RooCruijff(void *p);
   static void destruct_RooCruijff(void *p);
   static void streamer_RooCruijff(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooCruijff*)
   {
      ::RooCruijff *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooCruijff >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooCruijff", ::RooCruijff::Class_Version(), "RooCruijff.hh", 20,
                  typeid(::RooCruijff), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooCruijff::Dictionary, isa_proxy, 16,
                  sizeof(::RooCruijff) );
      instance.SetDelete(&delete_RooCruijff);
      instance.SetDeleteArray(&deleteArray_RooCruijff);
      instance.SetDestructor(&destruct_RooCruijff);
      instance.SetStreamerFunc(&streamer_RooCruijff);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooCruijff*)
   {
      return GenerateInitInstanceLocal((::RooCruijff*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RooCruijff*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarBinned(void *p = 0);
   static void *newArray_rarBinned(Long_t size, void *p);
   static void delete_rarBinned(void *p);
   static void deleteArray_rarBinned(void *p);
   static void destruct_rarBinned(void *p);
   static void streamer_rarBinned(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarBinned*)
   {
      ::rarBinned *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarBinned >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarBinned", ::rarBinned::Class_Version(), "rarBinned.hh", 40,
                  typeid(::rarBinned), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarBinned::Dictionary, isa_proxy, 16,
                  sizeof(::rarBinned) );
      instance.SetNew(&new_rarBinned);
      instance.SetNewArray(&newArray_rarBinned);
      instance.SetDelete(&delete_rarBinned);
      instance.SetDeleteArray(&deleteArray_rarBinned);
      instance.SetDestructor(&destruct_rarBinned);
      instance.SetStreamerFunc(&streamer_rarBinned);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarBinned*)
   {
      return GenerateInitInstanceLocal((::rarBinned*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarBinned*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarBallack(void *p = 0);
   static void *newArray_rarBallack(Long_t size, void *p);
   static void delete_rarBallack(void *p);
   static void deleteArray_rarBallack(void *p);
   static void destruct_rarBallack(void *p);
   static void streamer_rarBallack(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarBallack*)
   {
      ::rarBallack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarBallack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarBallack", ::rarBallack::Class_Version(), "rarBallack.hh", 40,
                  typeid(::rarBallack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarBallack::Dictionary, isa_proxy, 16,
                  sizeof(::rarBallack) );
      instance.SetNew(&new_rarBallack);
      instance.SetNewArray(&newArray_rarBallack);
      instance.SetDelete(&delete_rarBallack);
      instance.SetDeleteArray(&deleteArray_rarBallack);
      instance.SetDestructor(&destruct_rarBallack);
      instance.SetStreamerFunc(&streamer_rarBallack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarBallack*)
   {
      return GenerateInitInstanceLocal((::rarBallack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarBallack*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarRelBreitWigner(void *p = 0);
   static void *newArray_rarRelBreitWigner(Long_t size, void *p);
   static void delete_rarRelBreitWigner(void *p);
   static void deleteArray_rarRelBreitWigner(void *p);
   static void destruct_rarRelBreitWigner(void *p);
   static void streamer_rarRelBreitWigner(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarRelBreitWigner*)
   {
      ::rarRelBreitWigner *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarRelBreitWigner >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarRelBreitWigner", ::rarRelBreitWigner::Class_Version(), "rarRelBreitWigner.hh", 37,
                  typeid(::rarRelBreitWigner), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarRelBreitWigner::Dictionary, isa_proxy, 16,
                  sizeof(::rarRelBreitWigner) );
      instance.SetNew(&new_rarRelBreitWigner);
      instance.SetNewArray(&newArray_rarRelBreitWigner);
      instance.SetDelete(&delete_rarRelBreitWigner);
      instance.SetDeleteArray(&deleteArray_rarRelBreitWigner);
      instance.SetDestructor(&destruct_rarRelBreitWigner);
      instance.SetStreamerFunc(&streamer_rarRelBreitWigner);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarRelBreitWigner*)
   {
      return GenerateInitInstanceLocal((::rarRelBreitWigner*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarRelBreitWigner*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarBifurGauss(void *p = 0);
   static void *newArray_rarBifurGauss(Long_t size, void *p);
   static void delete_rarBifurGauss(void *p);
   static void deleteArray_rarBifurGauss(void *p);
   static void destruct_rarBifurGauss(void *p);
   static void streamer_rarBifurGauss(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarBifurGauss*)
   {
      ::rarBifurGauss *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarBifurGauss >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarBifurGauss", ::rarBifurGauss::Class_Version(), "rarBifurGauss.hh", 50,
                  typeid(::rarBifurGauss), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarBifurGauss::Dictionary, isa_proxy, 16,
                  sizeof(::rarBifurGauss) );
      instance.SetNew(&new_rarBifurGauss);
      instance.SetNewArray(&newArray_rarBifurGauss);
      instance.SetDelete(&delete_rarBifurGauss);
      instance.SetDeleteArray(&deleteArray_rarBifurGauss);
      instance.SetDestructor(&destruct_rarBifurGauss);
      instance.SetStreamerFunc(&streamer_rarBifurGauss);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarBifurGauss*)
   {
      return GenerateInitInstanceLocal((::rarBifurGauss*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarBifurGauss*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarUsrPdf(void *p = 0);
   static void *newArray_rarUsrPdf(Long_t size, void *p);
   static void delete_rarUsrPdf(void *p);
   static void deleteArray_rarUsrPdf(void *p);
   static void destruct_rarUsrPdf(void *p);
   static void streamer_rarUsrPdf(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarUsrPdf*)
   {
      ::rarUsrPdf *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarUsrPdf >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarUsrPdf", ::rarUsrPdf::Class_Version(), "rarUsrPdf.hh", 24,
                  typeid(::rarUsrPdf), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarUsrPdf::Dictionary, isa_proxy, 16,
                  sizeof(::rarUsrPdf) );
      instance.SetNew(&new_rarUsrPdf);
      instance.SetNewArray(&newArray_rarUsrPdf);
      instance.SetDelete(&delete_rarUsrPdf);
      instance.SetDeleteArray(&deleteArray_rarUsrPdf);
      instance.SetDestructor(&destruct_rarUsrPdf);
      instance.SetStreamerFunc(&streamer_rarUsrPdf);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarUsrPdf*)
   {
      return GenerateInitInstanceLocal((::rarUsrPdf*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarUsrPdf*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_RooGounarisSakurai(void *p);
   static void deleteArray_RooGounarisSakurai(void *p);
   static void destruct_RooGounarisSakurai(void *p);
   static void streamer_RooGounarisSakurai(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooGounarisSakurai*)
   {
      ::RooGounarisSakurai *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooGounarisSakurai >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooGounarisSakurai", ::RooGounarisSakurai::Class_Version(), "RooGounarisSakurai.hh", 16,
                  typeid(::RooGounarisSakurai), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooGounarisSakurai::Dictionary, isa_proxy, 16,
                  sizeof(::RooGounarisSakurai) );
      instance.SetDelete(&delete_RooGounarisSakurai);
      instance.SetDeleteArray(&deleteArray_RooGounarisSakurai);
      instance.SetDestructor(&destruct_RooGounarisSakurai);
      instance.SetStreamerFunc(&streamer_RooGounarisSakurai);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooGounarisSakurai*)
   {
      return GenerateInitInstanceLocal((::RooGounarisSakurai*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RooGounarisSakurai*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarCBShape(void *p = 0);
   static void *newArray_rarCBShape(Long_t size, void *p);
   static void delete_rarCBShape(void *p);
   static void deleteArray_rarCBShape(void *p);
   static void destruct_rarCBShape(void *p);
   static void streamer_rarCBShape(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarCBShape*)
   {
      ::rarCBShape *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarCBShape >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarCBShape", ::rarCBShape::Class_Version(), "rarCBShape.hh", 37,
                  typeid(::rarCBShape), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarCBShape::Dictionary, isa_proxy, 16,
                  sizeof(::rarCBShape) );
      instance.SetNew(&new_rarCBShape);
      instance.SetNewArray(&newArray_rarCBShape);
      instance.SetDelete(&delete_rarCBShape);
      instance.SetDeleteArray(&deleteArray_rarCBShape);
      instance.SetDestructor(&destruct_rarCBShape);
      instance.SetStreamerFunc(&streamer_rarCBShape);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarCBShape*)
   {
      return GenerateInitInstanceLocal((::rarCBShape*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarCBShape*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarGeneric(void *p = 0);
   static void *newArray_rarGeneric(Long_t size, void *p);
   static void delete_rarGeneric(void *p);
   static void deleteArray_rarGeneric(void *p);
   static void destruct_rarGeneric(void *p);
   static void streamer_rarGeneric(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarGeneric*)
   {
      ::rarGeneric *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarGeneric >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarGeneric", ::rarGeneric::Class_Version(), "rarGeneric.hh", 26,
                  typeid(::rarGeneric), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarGeneric::Dictionary, isa_proxy, 16,
                  sizeof(::rarGeneric) );
      instance.SetNew(&new_rarGeneric);
      instance.SetNewArray(&newArray_rarGeneric);
      instance.SetDelete(&delete_rarGeneric);
      instance.SetDeleteArray(&deleteArray_rarGeneric);
      instance.SetDestructor(&destruct_rarGeneric);
      instance.SetStreamerFunc(&streamer_rarGeneric);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarGeneric*)
   {
      return GenerateInitInstanceLocal((::rarGeneric*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarGeneric*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_RooBinnedPdf(void *p);
   static void deleteArray_RooBinnedPdf(void *p);
   static void destruct_RooBinnedPdf(void *p);
   static void streamer_RooBinnedPdf(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooBinnedPdf*)
   {
      ::RooBinnedPdf *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooBinnedPdf >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooBinnedPdf", ::RooBinnedPdf::Class_Version(), "RooBinnedPdf.hh", 30,
                  typeid(::RooBinnedPdf), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooBinnedPdf::Dictionary, isa_proxy, 16,
                  sizeof(::RooBinnedPdf) );
      instance.SetDelete(&delete_RooBinnedPdf);
      instance.SetDeleteArray(&deleteArray_RooBinnedPdf);
      instance.SetDestructor(&destruct_RooBinnedPdf);
      instance.SetStreamerFunc(&streamer_RooBinnedPdf);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooBinnedPdf*)
   {
      return GenerateInitInstanceLocal((::RooBinnedPdf*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RooBinnedPdf*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarCruijff(void *p = 0);
   static void *newArray_rarCruijff(Long_t size, void *p);
   static void delete_rarCruijff(void *p);
   static void deleteArray_rarCruijff(void *p);
   static void destruct_rarCruijff(void *p);
   static void streamer_rarCruijff(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarCruijff*)
   {
      ::rarCruijff *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarCruijff >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarCruijff", ::rarCruijff::Class_Version(), "rarCruijff.hh", 24,
                  typeid(::rarCruijff), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarCruijff::Dictionary, isa_proxy, 16,
                  sizeof(::rarCruijff) );
      instance.SetNew(&new_rarCruijff);
      instance.SetNewArray(&newArray_rarCruijff);
      instance.SetDelete(&delete_rarCruijff);
      instance.SetDeleteArray(&deleteArray_rarCruijff);
      instance.SetDestructor(&destruct_rarCruijff);
      instance.SetStreamerFunc(&streamer_rarCruijff);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarCruijff*)
   {
      return GenerateInitInstanceLocal((::rarCruijff*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarCruijff*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_rarMinuit(void *p);
   static void deleteArray_rarMinuit(void *p);
   static void destruct_rarMinuit(void *p);
   static void streamer_rarMinuit(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarMinuit*)
   {
      ::rarMinuit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarMinuit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarMinuit", ::rarMinuit::Class_Version(), "rarMinuit.hh", 27,
                  typeid(::rarMinuit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarMinuit::Dictionary, isa_proxy, 16,
                  sizeof(::rarMinuit) );
      instance.SetDelete(&delete_rarMinuit);
      instance.SetDeleteArray(&deleteArray_rarMinuit);
      instance.SetDestructor(&destruct_rarMinuit);
      instance.SetStreamerFunc(&streamer_rarMinuit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarMinuit*)
   {
      return GenerateInitInstanceLocal((::rarMinuit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarMinuit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarSimPdf(void *p = 0);
   static void *newArray_rarSimPdf(Long_t size, void *p);
   static void delete_rarSimPdf(void *p);
   static void deleteArray_rarSimPdf(void *p);
   static void destruct_rarSimPdf(void *p);
   static void streamer_rarSimPdf(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarSimPdf*)
   {
      ::rarSimPdf *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarSimPdf >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarSimPdf", ::rarSimPdf::Class_Version(), "rarSimPdf.hh", 29,
                  typeid(::rarSimPdf), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarSimPdf::Dictionary, isa_proxy, 16,
                  sizeof(::rarSimPdf) );
      instance.SetNew(&new_rarSimPdf);
      instance.SetNewArray(&newArray_rarSimPdf);
      instance.SetDelete(&delete_rarSimPdf);
      instance.SetDeleteArray(&deleteArray_rarSimPdf);
      instance.SetDestructor(&destruct_rarSimPdf);
      instance.SetStreamerFunc(&streamer_rarSimPdf);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarSimPdf*)
   {
      return GenerateInitInstanceLocal((::rarSimPdf*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarSimPdf*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarFlatte(void *p = 0);
   static void *newArray_rarFlatte(Long_t size, void *p);
   static void delete_rarFlatte(void *p);
   static void deleteArray_rarFlatte(void *p);
   static void destruct_rarFlatte(void *p);
   static void streamer_rarFlatte(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarFlatte*)
   {
      ::rarFlatte *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarFlatte >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarFlatte", ::rarFlatte::Class_Version(), "rarFlatte.hh", 46,
                  typeid(::rarFlatte), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarFlatte::Dictionary, isa_proxy, 16,
                  sizeof(::rarFlatte) );
      instance.SetNew(&new_rarFlatte);
      instance.SetNewArray(&newArray_rarFlatte);
      instance.SetDelete(&delete_rarFlatte);
      instance.SetDeleteArray(&deleteArray_rarFlatte);
      instance.SetDestructor(&destruct_rarFlatte);
      instance.SetStreamerFunc(&streamer_rarFlatte);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarFlatte*)
   {
      return GenerateInitInstanceLocal((::rarFlatte*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarFlatte*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarVoigtian(void *p = 0);
   static void *newArray_rarVoigtian(Long_t size, void *p);
   static void delete_rarVoigtian(void *p);
   static void deleteArray_rarVoigtian(void *p);
   static void destruct_rarVoigtian(void *p);
   static void streamer_rarVoigtian(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarVoigtian*)
   {
      ::rarVoigtian *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarVoigtian >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarVoigtian", ::rarVoigtian::Class_Version(), "rarVoigtian.hh", 35,
                  typeid(::rarVoigtian), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarVoigtian::Dictionary, isa_proxy, 16,
                  sizeof(::rarVoigtian) );
      instance.SetNew(&new_rarVoigtian);
      instance.SetNewArray(&newArray_rarVoigtian);
      instance.SetDelete(&delete_rarVoigtian);
      instance.SetDeleteArray(&deleteArray_rarVoigtian);
      instance.SetDestructor(&destruct_rarVoigtian);
      instance.SetStreamerFunc(&streamer_rarVoigtian);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarVoigtian*)
   {
      return GenerateInitInstanceLocal((::rarVoigtian*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarVoigtian*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarThreshold(void *p = 0);
   static void *newArray_rarThreshold(Long_t size, void *p);
   static void delete_rarThreshold(void *p);
   static void deleteArray_rarThreshold(void *p);
   static void destruct_rarThreshold(void *p);
   static void streamer_rarThreshold(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarThreshold*)
   {
      ::rarThreshold *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarThreshold >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarThreshold", ::rarThreshold::Class_Version(), "rarThreshold.hh", 41,
                  typeid(::rarThreshold), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarThreshold::Dictionary, isa_proxy, 16,
                  sizeof(::rarThreshold) );
      instance.SetNew(&new_rarThreshold);
      instance.SetNewArray(&newArray_rarThreshold);
      instance.SetDelete(&delete_rarThreshold);
      instance.SetDeleteArray(&deleteArray_rarThreshold);
      instance.SetDestructor(&destruct_rarThreshold);
      instance.SetStreamerFunc(&streamer_rarThreshold);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarThreshold*)
   {
      return GenerateInitInstanceLocal((::rarThreshold*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarThreshold*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarGaussModel(void *p = 0);
   static void *newArray_rarGaussModel(Long_t size, void *p);
   static void delete_rarGaussModel(void *p);
   static void deleteArray_rarGaussModel(void *p);
   static void destruct_rarGaussModel(void *p);
   static void streamer_rarGaussModel(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarGaussModel*)
   {
      ::rarGaussModel *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarGaussModel >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarGaussModel", ::rarGaussModel::Class_Version(), "rarGaussModel.hh", 48,
                  typeid(::rarGaussModel), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarGaussModel::Dictionary, isa_proxy, 16,
                  sizeof(::rarGaussModel) );
      instance.SetNew(&new_rarGaussModel);
      instance.SetNewArray(&newArray_rarGaussModel);
      instance.SetDelete(&delete_rarGaussModel);
      instance.SetDeleteArray(&deleteArray_rarGaussModel);
      instance.SetDestructor(&destruct_rarGaussModel);
      instance.SetStreamerFunc(&streamer_rarGaussModel);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarGaussModel*)
   {
      return GenerateInitInstanceLocal((::rarGaussModel*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarGaussModel*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarExp(void *p = 0);
   static void *newArray_rarExp(Long_t size, void *p);
   static void delete_rarExp(void *p);
   static void deleteArray_rarExp(void *p);
   static void destruct_rarExp(void *p);
   static void streamer_rarExp(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarExp*)
   {
      ::rarExp *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarExp >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarExp", ::rarExp::Class_Version(), "rarExp.hh", 31,
                  typeid(::rarExp), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarExp::Dictionary, isa_proxy, 16,
                  sizeof(::rarExp) );
      instance.SetNew(&new_rarExp);
      instance.SetNewArray(&newArray_rarExp);
      instance.SetDelete(&delete_rarExp);
      instance.SetDeleteArray(&deleteArray_rarExp);
      instance.SetDestructor(&destruct_rarExp);
      instance.SetStreamerFunc(&streamer_rarExp);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarExp*)
   {
      return GenerateInitInstanceLocal((::rarExp*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarExp*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarToyList(void *p = 0);
   static void *newArray_rarToyList(Long_t size, void *p);
   static void delete_rarToyList(void *p);
   static void deleteArray_rarToyList(void *p);
   static void destruct_rarToyList(void *p);
   static void streamer_rarToyList(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarToyList*)
   {
      ::rarToyList *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarToyList >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarToyList", ::rarToyList::Class_Version(), "rarToyList.hh", 25,
                  typeid(::rarToyList), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarToyList::Dictionary, isa_proxy, 16,
                  sizeof(::rarToyList) );
      instance.SetNew(&new_rarToyList);
      instance.SetNewArray(&newArray_rarToyList);
      instance.SetDelete(&delete_rarToyList);
      instance.SetDeleteArray(&deleteArray_rarToyList);
      instance.SetDestructor(&destruct_rarToyList);
      instance.SetStreamerFunc(&streamer_rarToyList);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarToyList*)
   {
      return GenerateInitInstanceLocal((::rarToyList*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarToyList*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarArgusBG(void *p = 0);
   static void *newArray_rarArgusBG(Long_t size, void *p);
   static void delete_rarArgusBG(void *p);
   static void deleteArray_rarArgusBG(void *p);
   static void destruct_rarArgusBG(void *p);
   static void streamer_rarArgusBG(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarArgusBG*)
   {
      ::rarArgusBG *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarArgusBG >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarArgusBG", ::rarArgusBG::Class_Version(), "rarArgusBG.hh", 36,
                  typeid(::rarArgusBG), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarArgusBG::Dictionary, isa_proxy, 16,
                  sizeof(::rarArgusBG) );
      instance.SetNew(&new_rarArgusBG);
      instance.SetNewArray(&newArray_rarArgusBG);
      instance.SetDelete(&delete_rarArgusBG);
      instance.SetDeleteArray(&deleteArray_rarArgusBG);
      instance.SetDestructor(&destruct_rarArgusBG);
      instance.SetStreamerFunc(&streamer_rarArgusBG);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarArgusBG*)
   {
      return GenerateInitInstanceLocal((::rarArgusBG*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarArgusBG*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_RooRelBreitWigner(void *p);
   static void deleteArray_RooRelBreitWigner(void *p);
   static void destruct_RooRelBreitWigner(void *p);
   static void streamer_RooRelBreitWigner(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooRelBreitWigner*)
   {
      ::RooRelBreitWigner *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooRelBreitWigner >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooRelBreitWigner", ::RooRelBreitWigner::Class_Version(), "RooRelBreitWigner.hh", 9,
                  typeid(::RooRelBreitWigner), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooRelBreitWigner::Dictionary, isa_proxy, 16,
                  sizeof(::RooRelBreitWigner) );
      instance.SetDelete(&delete_RooRelBreitWigner);
      instance.SetDeleteArray(&deleteArray_RooRelBreitWigner);
      instance.SetDestructor(&destruct_RooRelBreitWigner);
      instance.SetStreamerFunc(&streamer_RooRelBreitWigner);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooRelBreitWigner*)
   {
      return GenerateInitInstanceLocal((::RooRelBreitWigner*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RooRelBreitWigner*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarGounarisSakurai(void *p = 0);
   static void *newArray_rarGounarisSakurai(Long_t size, void *p);
   static void delete_rarGounarisSakurai(void *p);
   static void deleteArray_rarGounarisSakurai(void *p);
   static void destruct_rarGounarisSakurai(void *p);
   static void streamer_rarGounarisSakurai(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarGounarisSakurai*)
   {
      ::rarGounarisSakurai *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarGounarisSakurai >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarGounarisSakurai", ::rarGounarisSakurai::Class_Version(), "rarGounarisSakurai.hh", 45,
                  typeid(::rarGounarisSakurai), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarGounarisSakurai::Dictionary, isa_proxy, 16,
                  sizeof(::rarGounarisSakurai) );
      instance.SetNew(&new_rarGounarisSakurai);
      instance.SetNewArray(&newArray_rarGounarisSakurai);
      instance.SetDelete(&delete_rarGounarisSakurai);
      instance.SetDeleteArray(&deleteArray_rarGounarisSakurai);
      instance.SetDestructor(&destruct_rarGounarisSakurai);
      instance.SetStreamerFunc(&streamer_rarGounarisSakurai);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarGounarisSakurai*)
   {
      return GenerateInitInstanceLocal((::rarGounarisSakurai*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarGounarisSakurai*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarMLFitter(void *p = 0);
   static void *newArray_rarMLFitter(Long_t size, void *p);
   static void delete_rarMLFitter(void *p);
   static void deleteArray_rarMLFitter(void *p);
   static void destruct_rarMLFitter(void *p);
   static void streamer_rarMLFitter(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarMLFitter*)
   {
      ::rarMLFitter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarMLFitter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarMLFitter", ::rarMLFitter::Class_Version(), "rarMLFitter.hh", 49,
                  typeid(::rarMLFitter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarMLFitter::Dictionary, isa_proxy, 16,
                  sizeof(::rarMLFitter) );
      instance.SetNew(&new_rarMLFitter);
      instance.SetNewArray(&newArray_rarMLFitter);
      instance.SetDelete(&delete_rarMLFitter);
      instance.SetDeleteArray(&deleteArray_rarMLFitter);
      instance.SetDestructor(&destruct_rarMLFitter);
      instance.SetStreamerFunc(&streamer_rarMLFitter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarMLFitter*)
   {
      return GenerateInitInstanceLocal((::rarMLFitter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarMLFitter*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarKeys(void *p = 0);
   static void *newArray_rarKeys(Long_t size, void *p);
   static void delete_rarKeys(void *p);
   static void deleteArray_rarKeys(void *p);
   static void destruct_rarKeys(void *p);
   static void streamer_rarKeys(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarKeys*)
   {
      ::rarKeys *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarKeys >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarKeys", ::rarKeys::Class_Version(), "rarKeys.hh", 45,
                  typeid(::rarKeys), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarKeys::Dictionary, isa_proxy, 16,
                  sizeof(::rarKeys) );
      instance.SetNew(&new_rarKeys);
      instance.SetNewArray(&newArray_rarKeys);
      instance.SetDelete(&delete_rarKeys);
      instance.SetDeleteArray(&deleteArray_rarKeys);
      instance.SetDestructor(&destruct_rarKeys);
      instance.SetStreamerFunc(&streamer_rarKeys);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarKeys*)
   {
      return GenerateInitInstanceLocal((::rarKeys*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarKeys*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarMLPdf(void *p = 0);
   static void *newArray_rarMLPdf(Long_t size, void *p);
   static void delete_rarMLPdf(void *p);
   static void deleteArray_rarMLPdf(void *p);
   static void destruct_rarMLPdf(void *p);
   static void streamer_rarMLPdf(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarMLPdf*)
   {
      ::rarMLPdf *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarMLPdf >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarMLPdf", ::rarMLPdf::Class_Version(), "rarMLPdf.hh", 43,
                  typeid(::rarMLPdf), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarMLPdf::Dictionary, isa_proxy, 16,
                  sizeof(::rarMLPdf) );
      instance.SetNew(&new_rarMLPdf);
      instance.SetNewArray(&newArray_rarMLPdf);
      instance.SetDelete(&delete_rarMLPdf);
      instance.SetDeleteArray(&deleteArray_rarMLPdf);
      instance.SetDestructor(&destruct_rarMLPdf);
      instance.SetStreamerFunc(&streamer_rarMLPdf);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarMLPdf*)
   {
      return GenerateInitInstanceLocal((::rarMLPdf*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarMLPdf*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_RooFlatte(void *p);
   static void deleteArray_RooFlatte(void *p);
   static void destruct_RooFlatte(void *p);
   static void streamer_RooFlatte(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooFlatte*)
   {
      ::RooFlatte *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooFlatte >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooFlatte", ::RooFlatte::Class_Version(), "RooFlatte.hh", 9,
                  typeid(::RooFlatte), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooFlatte::Dictionary, isa_proxy, 16,
                  sizeof(::RooFlatte) );
      instance.SetDelete(&delete_RooFlatte);
      instance.SetDeleteArray(&deleteArray_RooFlatte);
      instance.SetDestructor(&destruct_RooFlatte);
      instance.SetStreamerFunc(&streamer_RooFlatte);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooFlatte*)
   {
      return GenerateInitInstanceLocal((::RooFlatte*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RooFlatte*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarDecay(void *p = 0);
   static void *newArray_rarDecay(Long_t size, void *p);
   static void delete_rarDecay(void *p);
   static void deleteArray_rarDecay(void *p);
   static void destruct_rarDecay(void *p);
   static void streamer_rarDecay(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarDecay*)
   {
      ::rarDecay *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarDecay >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarDecay", ::rarDecay::Class_Version(), "rarDecay.hh", 31,
                  typeid(::rarDecay), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarDecay::Dictionary, isa_proxy, 16,
                  sizeof(::rarDecay) );
      instance.SetNew(&new_rarDecay);
      instance.SetNewArray(&newArray_rarDecay);
      instance.SetDelete(&delete_rarDecay);
      instance.SetDeleteArray(&deleteArray_rarDecay);
      instance.SetDestructor(&destruct_rarDecay);
      instance.SetStreamerFunc(&streamer_rarDecay);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarDecay*)
   {
      return GenerateInitInstanceLocal((::rarDecay*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarDecay*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_rarPoly(void *p = 0);
   static void *newArray_rarPoly(Long_t size, void *p);
   static void delete_rarPoly(void *p);
   static void deleteArray_rarPoly(void *p);
   static void destruct_rarPoly(void *p);
   static void streamer_rarPoly(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::rarPoly*)
   {
      ::rarPoly *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::rarPoly >(0);
      static ::ROOT::TGenericClassInfo 
         instance("rarPoly", ::rarPoly::Class_Version(), "rarPoly.hh", 41,
                  typeid(::rarPoly), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::rarPoly::Dictionary, isa_proxy, 16,
                  sizeof(::rarPoly) );
      instance.SetNew(&new_rarPoly);
      instance.SetNewArray(&newArray_rarPoly);
      instance.SetDelete(&delete_rarPoly);
      instance.SetDeleteArray(&deleteArray_rarPoly);
      instance.SetDestructor(&destruct_rarPoly);
      instance.SetStreamerFunc(&streamer_rarPoly);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::rarPoly*)
   {
      return GenerateInitInstanceLocal((::rarPoly*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::rarPoly*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr rarStrParser::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarStrParser::Class_Name()
{
   return "rarStrParser";
}

//______________________________________________________________________________
const char *rarStrParser::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarStrParser*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarStrParser::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarStrParser*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarStrParser::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarStrParser*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarStrParser::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarStrParser*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RooLass::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooLass::Class_Name()
{
   return "RooLass";
}

//______________________________________________________________________________
const char *RooLass::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooLass*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooLass::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooLass*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooLass::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooLass*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooLass::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooLass*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarConfig::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarConfig::Class_Name()
{
   return "rarConfig";
}

//______________________________________________________________________________
const char *rarConfig::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarConfig*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarConfig::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarConfig*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarConfig::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarConfig*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarConfig::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarConfig*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarDatasetDef::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarDatasetDef::Class_Name()
{
   return "rarDatasetDef";
}

//______________________________________________________________________________
const char *rarDatasetDef::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarDatasetDef*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarDatasetDef::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarDatasetDef*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarDatasetDef::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarDatasetDef*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarDatasetDef::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarDatasetDef*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarDatasets::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarDatasets::Class_Name()
{
   return "rarDatasets";
}

//______________________________________________________________________________
const char *rarDatasets::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarDatasets*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarDatasets::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarDatasets*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarDatasets::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarDatasets*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarDatasets::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarDatasets*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarBasePdf::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarBasePdf::Class_Name()
{
   return "rarBasePdf";
}

//______________________________________________________________________________
const char *rarBasePdf::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarBasePdf*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarBasePdf::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarBasePdf*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarBasePdf::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarBasePdf*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarBasePdf::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarBasePdf*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarTwoGauss::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarTwoGauss::Class_Name()
{
   return "rarTwoGauss";
}

//______________________________________________________________________________
const char *rarTwoGauss::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarTwoGauss*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarTwoGauss::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarTwoGauss*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarTwoGauss::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarTwoGauss*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarTwoGauss::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarTwoGauss*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarCompBase::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarCompBase::Class_Name()
{
   return "rarCompBase";
}

//______________________________________________________________________________
const char *rarCompBase::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarCompBase*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarCompBase::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarCompBase*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarCompBase::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarCompBase*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarCompBase::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarCompBase*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarMultPdf::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarMultPdf::Class_Name()
{
   return "rarMultPdf";
}

//______________________________________________________________________________
const char *rarMultPdf::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarMultPdf*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarMultPdf::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarMultPdf*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarMultPdf::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarMultPdf*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarMultPdf::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarMultPdf*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RooOsipDisc::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooOsipDisc::Class_Name()
{
   return "RooOsipDisc";
}

//______________________________________________________________________________
const char *RooOsipDisc::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooOsipDisc*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooOsipDisc::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooOsipDisc*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooOsipDisc::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooOsipDisc*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooOsipDisc::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooOsipDisc*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarStep::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarStep::Class_Name()
{
   return "rarStep";
}

//______________________________________________________________________________
const char *rarStep::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarStep*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarStep::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarStep*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarStep::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarStep*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarStep::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarStep*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarProd::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarProd::Class_Name()
{
   return "rarProd";
}

//______________________________________________________________________________
const char *rarProd::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarProd*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarProd::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarProd*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarProd::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarProd*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarProd::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarProd*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarNLL::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarNLL::Class_Name()
{
   return "rarNLL";
}

//______________________________________________________________________________
const char *rarNLL::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarNLL*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarNLL::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarNLL*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarNLL::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarNLL*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarNLL::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarNLL*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarUniform::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarUniform::Class_Name()
{
   return "rarUniform";
}

//______________________________________________________________________________
const char *rarUniform::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarUniform*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarUniform::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarUniform*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarUniform::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarUniform*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarUniform::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarUniform*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarOsipDisc::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarOsipDisc::Class_Name()
{
   return "rarOsipDisc";
}

//______________________________________________________________________________
const char *rarOsipDisc::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarOsipDisc*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarOsipDisc::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarOsipDisc*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarOsipDisc::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarOsipDisc*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarOsipDisc::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarOsipDisc*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RooBallack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooBallack::Class_Name()
{
   return "RooBallack";
}

//______________________________________________________________________________
const char *RooBallack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooBallack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooBallack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooBallack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooBallack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooBallack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooBallack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooBallack*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarLass::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarLass::Class_Name()
{
   return "rarLass";
}

//______________________________________________________________________________
const char *rarLass::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarLass*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarLass::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarLass*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarLass::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarLass*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarLass::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarLass*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarAdd::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarAdd::Class_Name()
{
   return "rarAdd";
}

//______________________________________________________________________________
const char *rarAdd::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarAdd*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarAdd::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarAdd*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarAdd::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarAdd*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarAdd::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarAdd*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarTriGauss::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarTriGauss::Class_Name()
{
   return "rarTriGauss";
}

//______________________________________________________________________________
const char *rarTriGauss::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarTriGauss*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarTriGauss::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarTriGauss*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarTriGauss::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarTriGauss*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarTriGauss::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarTriGauss*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarNovosibirsk::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarNovosibirsk::Class_Name()
{
   return "rarNovosibirsk";
}

//______________________________________________________________________________
const char *rarNovosibirsk::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarNovosibirsk*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarNovosibirsk::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarNovosibirsk*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarNovosibirsk::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarNovosibirsk*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarNovosibirsk::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarNovosibirsk*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarHistPdf::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarHistPdf::Class_Name()
{
   return "rarHistPdf";
}

//______________________________________________________________________________
const char *rarHistPdf::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarHistPdf*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarHistPdf::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarHistPdf*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarHistPdf::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarHistPdf*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarHistPdf::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarHistPdf*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarSPlot::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarSPlot::Class_Name()
{
   return "rarSPlot";
}

//______________________________________________________________________________
const char *rarSPlot::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarSPlot*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarSPlot::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarSPlot*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarSPlot::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarSPlot*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarSPlot::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarSPlot*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RooThreshold::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooThreshold::Class_Name()
{
   return "RooThreshold";
}

//______________________________________________________________________________
const char *RooThreshold::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooThreshold*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooThreshold::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooThreshold*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooThreshold::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooThreshold*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooThreshold::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooThreshold*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarGaussian::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarGaussian::Class_Name()
{
   return "rarGaussian";
}

//______________________________________________________________________________
const char *rarGaussian::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarGaussian*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarGaussian::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarGaussian*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarGaussian::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarGaussian*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarGaussian::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarGaussian*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RooCruijff::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooCruijff::Class_Name()
{
   return "RooCruijff";
}

//______________________________________________________________________________
const char *RooCruijff::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooCruijff*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooCruijff::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooCruijff*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooCruijff::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooCruijff*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooCruijff::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooCruijff*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarBinned::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarBinned::Class_Name()
{
   return "rarBinned";
}

//______________________________________________________________________________
const char *rarBinned::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarBinned*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarBinned::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarBinned*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarBinned::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarBinned*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarBinned::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarBinned*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarBallack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarBallack::Class_Name()
{
   return "rarBallack";
}

//______________________________________________________________________________
const char *rarBallack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarBallack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarBallack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarBallack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarBallack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarBallack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarBallack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarBallack*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarRelBreitWigner::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarRelBreitWigner::Class_Name()
{
   return "rarRelBreitWigner";
}

//______________________________________________________________________________
const char *rarRelBreitWigner::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarRelBreitWigner*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarRelBreitWigner::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarRelBreitWigner*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarRelBreitWigner::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarRelBreitWigner*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarRelBreitWigner::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarRelBreitWigner*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarBifurGauss::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarBifurGauss::Class_Name()
{
   return "rarBifurGauss";
}

//______________________________________________________________________________
const char *rarBifurGauss::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarBifurGauss*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarBifurGauss::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarBifurGauss*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarBifurGauss::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarBifurGauss*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarBifurGauss::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarBifurGauss*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarUsrPdf::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarUsrPdf::Class_Name()
{
   return "rarUsrPdf";
}

//______________________________________________________________________________
const char *rarUsrPdf::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarUsrPdf*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarUsrPdf::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarUsrPdf*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarUsrPdf::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarUsrPdf*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarUsrPdf::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarUsrPdf*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RooGounarisSakurai::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooGounarisSakurai::Class_Name()
{
   return "RooGounarisSakurai";
}

//______________________________________________________________________________
const char *RooGounarisSakurai::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooGounarisSakurai*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooGounarisSakurai::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooGounarisSakurai*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooGounarisSakurai::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooGounarisSakurai*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooGounarisSakurai::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooGounarisSakurai*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarCBShape::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarCBShape::Class_Name()
{
   return "rarCBShape";
}

//______________________________________________________________________________
const char *rarCBShape::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarCBShape*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarCBShape::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarCBShape*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarCBShape::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarCBShape*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarCBShape::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarCBShape*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarGeneric::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarGeneric::Class_Name()
{
   return "rarGeneric";
}

//______________________________________________________________________________
const char *rarGeneric::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarGeneric*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarGeneric::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarGeneric*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarGeneric::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarGeneric*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarGeneric::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarGeneric*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RooBinnedPdf::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooBinnedPdf::Class_Name()
{
   return "RooBinnedPdf";
}

//______________________________________________________________________________
const char *RooBinnedPdf::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooBinnedPdf*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooBinnedPdf::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooBinnedPdf*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooBinnedPdf::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooBinnedPdf*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooBinnedPdf::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooBinnedPdf*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarCruijff::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarCruijff::Class_Name()
{
   return "rarCruijff";
}

//______________________________________________________________________________
const char *rarCruijff::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarCruijff*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarCruijff::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarCruijff*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarCruijff::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarCruijff*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarCruijff::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarCruijff*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarMinuit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarMinuit::Class_Name()
{
   return "rarMinuit";
}

//______________________________________________________________________________
const char *rarMinuit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarMinuit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarMinuit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarMinuit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarMinuit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarMinuit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarMinuit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarMinuit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarSimPdf::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarSimPdf::Class_Name()
{
   return "rarSimPdf";
}

//______________________________________________________________________________
const char *rarSimPdf::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarSimPdf*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarSimPdf::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarSimPdf*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarSimPdf::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarSimPdf*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarSimPdf::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarSimPdf*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarFlatte::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarFlatte::Class_Name()
{
   return "rarFlatte";
}

//______________________________________________________________________________
const char *rarFlatte::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarFlatte*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarFlatte::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarFlatte*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarFlatte::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarFlatte*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarFlatte::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarFlatte*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarVoigtian::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarVoigtian::Class_Name()
{
   return "rarVoigtian";
}

//______________________________________________________________________________
const char *rarVoigtian::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarVoigtian*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarVoigtian::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarVoigtian*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarVoigtian::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarVoigtian*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarVoigtian::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarVoigtian*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarThreshold::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarThreshold::Class_Name()
{
   return "rarThreshold";
}

//______________________________________________________________________________
const char *rarThreshold::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarThreshold*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarThreshold::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarThreshold*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarThreshold::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarThreshold*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarThreshold::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarThreshold*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarGaussModel::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarGaussModel::Class_Name()
{
   return "rarGaussModel";
}

//______________________________________________________________________________
const char *rarGaussModel::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarGaussModel*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarGaussModel::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarGaussModel*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarGaussModel::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarGaussModel*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarGaussModel::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarGaussModel*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarExp::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarExp::Class_Name()
{
   return "rarExp";
}

//______________________________________________________________________________
const char *rarExp::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarExp*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarExp::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarExp*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarExp::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarExp*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarExp::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarExp*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarToyList::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarToyList::Class_Name()
{
   return "rarToyList";
}

//______________________________________________________________________________
const char *rarToyList::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarToyList*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarToyList::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarToyList*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarToyList::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarToyList*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarToyList::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarToyList*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarArgusBG::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarArgusBG::Class_Name()
{
   return "rarArgusBG";
}

//______________________________________________________________________________
const char *rarArgusBG::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarArgusBG*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarArgusBG::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarArgusBG*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarArgusBG::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarArgusBG*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarArgusBG::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarArgusBG*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RooRelBreitWigner::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooRelBreitWigner::Class_Name()
{
   return "RooRelBreitWigner";
}

//______________________________________________________________________________
const char *RooRelBreitWigner::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooRelBreitWigner*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooRelBreitWigner::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooRelBreitWigner*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooRelBreitWigner::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooRelBreitWigner*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooRelBreitWigner::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooRelBreitWigner*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarGounarisSakurai::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarGounarisSakurai::Class_Name()
{
   return "rarGounarisSakurai";
}

//______________________________________________________________________________
const char *rarGounarisSakurai::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarGounarisSakurai*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarGounarisSakurai::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarGounarisSakurai*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarGounarisSakurai::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarGounarisSakurai*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarGounarisSakurai::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarGounarisSakurai*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarMLFitter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarMLFitter::Class_Name()
{
   return "rarMLFitter";
}

//______________________________________________________________________________
const char *rarMLFitter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarMLFitter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarMLFitter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarMLFitter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarMLFitter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarMLFitter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarMLFitter::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarMLFitter*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarKeys::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarKeys::Class_Name()
{
   return "rarKeys";
}

//______________________________________________________________________________
const char *rarKeys::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarKeys*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarKeys::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarKeys*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarKeys::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarKeys*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarKeys::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarKeys*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarMLPdf::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarMLPdf::Class_Name()
{
   return "rarMLPdf";
}

//______________________________________________________________________________
const char *rarMLPdf::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarMLPdf*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarMLPdf::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarMLPdf*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarMLPdf::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarMLPdf*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarMLPdf::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarMLPdf*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RooFlatte::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooFlatte::Class_Name()
{
   return "RooFlatte";
}

//______________________________________________________________________________
const char *RooFlatte::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooFlatte*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooFlatte::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooFlatte*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooFlatte::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooFlatte*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooFlatte::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooFlatte*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarDecay::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarDecay::Class_Name()
{
   return "rarDecay";
}

//______________________________________________________________________________
const char *rarDecay::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarDecay*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarDecay::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarDecay*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarDecay::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarDecay*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarDecay::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarDecay*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr rarPoly::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *rarPoly::Class_Name()
{
   return "rarPoly";
}

//______________________________________________________________________________
const char *rarPoly::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarPoly*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int rarPoly::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::rarPoly*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *rarPoly::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarPoly*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *rarPoly::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::rarPoly*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void rarStrParser::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarStrParser.

   TObject::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarStrParser(void *p) {
      return  p ? new(p) ::rarStrParser : new ::rarStrParser;
   }
   static void *newArray_rarStrParser(Long_t nElements, void *p) {
      return p ? new(p) ::rarStrParser[nElements] : new ::rarStrParser[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarStrParser(void *p) {
      delete ((::rarStrParser*)p);
   }
   static void deleteArray_rarStrParser(void *p) {
      delete [] ((::rarStrParser*)p);
   }
   static void destruct_rarStrParser(void *p) {
      typedef ::rarStrParser current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarStrParser(TBuffer &buf, void *obj) {
      ((::rarStrParser*)obj)->::rarStrParser::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarStrParser

//______________________________________________________________________________
void RooLass::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooLass.

   RooAbsPdf::Streamer(R__b);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RooLass(void *p) {
      delete ((::RooLass*)p);
   }
   static void deleteArray_RooLass(void *p) {
      delete [] ((::RooLass*)p);
   }
   static void destruct_RooLass(void *p) {
      typedef ::RooLass current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooLass(TBuffer &buf, void *obj) {
      ((::RooLass*)obj)->::RooLass::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooLass

//______________________________________________________________________________
void rarConfig::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarConfig.

   TNamed::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarConfig(void *p) {
      return  p ? new(p) ::rarConfig : new ::rarConfig;
   }
   static void *newArray_rarConfig(Long_t nElements, void *p) {
      return p ? new(p) ::rarConfig[nElements] : new ::rarConfig[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarConfig(void *p) {
      delete ((::rarConfig*)p);
   }
   static void deleteArray_rarConfig(void *p) {
      delete [] ((::rarConfig*)p);
   }
   static void destruct_rarConfig(void *p) {
      typedef ::rarConfig current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarConfig(TBuffer &buf, void *obj) {
      ((::rarConfig*)obj)->::rarConfig::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarConfig

//______________________________________________________________________________
void rarDatasetDef::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarDatasetDef.

   rarConfig::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarDatasetDef(void *p) {
      return  p ? new(p) ::rarDatasetDef : new ::rarDatasetDef;
   }
   static void *newArray_rarDatasetDef(Long_t nElements, void *p) {
      return p ? new(p) ::rarDatasetDef[nElements] : new ::rarDatasetDef[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarDatasetDef(void *p) {
      delete ((::rarDatasetDef*)p);
   }
   static void deleteArray_rarDatasetDef(void *p) {
      delete [] ((::rarDatasetDef*)p);
   }
   static void destruct_rarDatasetDef(void *p) {
      typedef ::rarDatasetDef current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarDatasetDef(TBuffer &buf, void *obj) {
      ((::rarDatasetDef*)obj)->::rarDatasetDef::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarDatasetDef

//______________________________________________________________________________
void rarDatasets::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarDatasets.

   rarConfig::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarDatasets(void *p) {
      return  p ? new(p) ::rarDatasets : new ::rarDatasets;
   }
   static void *newArray_rarDatasets(Long_t nElements, void *p) {
      return p ? new(p) ::rarDatasets[nElements] : new ::rarDatasets[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarDatasets(void *p) {
      delete ((::rarDatasets*)p);
   }
   static void deleteArray_rarDatasets(void *p) {
      delete [] ((::rarDatasets*)p);
   }
   static void destruct_rarDatasets(void *p) {
      typedef ::rarDatasets current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarDatasets(TBuffer &buf, void *obj) {
      ((::rarDatasets*)obj)->::rarDatasets::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarDatasets

//______________________________________________________________________________
void rarBasePdf::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarBasePdf.

   rarConfig::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarBasePdf(void *p) {
      return  p ? new(p) ::rarBasePdf : new ::rarBasePdf;
   }
   static void *newArray_rarBasePdf(Long_t nElements, void *p) {
      return p ? new(p) ::rarBasePdf[nElements] : new ::rarBasePdf[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarBasePdf(void *p) {
      delete ((::rarBasePdf*)p);
   }
   static void deleteArray_rarBasePdf(void *p) {
      delete [] ((::rarBasePdf*)p);
   }
   static void destruct_rarBasePdf(void *p) {
      typedef ::rarBasePdf current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarBasePdf(TBuffer &buf, void *obj) {
      ((::rarBasePdf*)obj)->::rarBasePdf::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarBasePdf

//______________________________________________________________________________
void rarTwoGauss::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarTwoGauss.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarTwoGauss(void *p) {
      return  p ? new(p) ::rarTwoGauss : new ::rarTwoGauss;
   }
   static void *newArray_rarTwoGauss(Long_t nElements, void *p) {
      return p ? new(p) ::rarTwoGauss[nElements] : new ::rarTwoGauss[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarTwoGauss(void *p) {
      delete ((::rarTwoGauss*)p);
   }
   static void deleteArray_rarTwoGauss(void *p) {
      delete [] ((::rarTwoGauss*)p);
   }
   static void destruct_rarTwoGauss(void *p) {
      typedef ::rarTwoGauss current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarTwoGauss(TBuffer &buf, void *obj) {
      ((::rarTwoGauss*)obj)->::rarTwoGauss::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarTwoGauss

//______________________________________________________________________________
void rarCompBase::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarCompBase.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarCompBase(void *p) {
      return  p ? new(p) ::rarCompBase : new ::rarCompBase;
   }
   static void *newArray_rarCompBase(Long_t nElements, void *p) {
      return p ? new(p) ::rarCompBase[nElements] : new ::rarCompBase[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarCompBase(void *p) {
      delete ((::rarCompBase*)p);
   }
   static void deleteArray_rarCompBase(void *p) {
      delete [] ((::rarCompBase*)p);
   }
   static void destruct_rarCompBase(void *p) {
      typedef ::rarCompBase current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarCompBase(TBuffer &buf, void *obj) {
      ((::rarCompBase*)obj)->::rarCompBase::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarCompBase

//______________________________________________________________________________
void rarMultPdf::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarMultPdf.

   rarCompBase::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarMultPdf(void *p) {
      return  p ? new(p) ::rarMultPdf : new ::rarMultPdf;
   }
   static void *newArray_rarMultPdf(Long_t nElements, void *p) {
      return p ? new(p) ::rarMultPdf[nElements] : new ::rarMultPdf[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarMultPdf(void *p) {
      delete ((::rarMultPdf*)p);
   }
   static void deleteArray_rarMultPdf(void *p) {
      delete [] ((::rarMultPdf*)p);
   }
   static void destruct_rarMultPdf(void *p) {
      typedef ::rarMultPdf current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarMultPdf(TBuffer &buf, void *obj) {
      ((::rarMultPdf*)obj)->::rarMultPdf::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarMultPdf

//______________________________________________________________________________
void RooOsipDisc::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooOsipDisc.

   RooAbsPdf::Streamer(R__b);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RooOsipDisc(void *p) {
      delete ((::RooOsipDisc*)p);
   }
   static void deleteArray_RooOsipDisc(void *p) {
      delete [] ((::RooOsipDisc*)p);
   }
   static void destruct_RooOsipDisc(void *p) {
      typedef ::RooOsipDisc current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooOsipDisc(TBuffer &buf, void *obj) {
      ((::RooOsipDisc*)obj)->::RooOsipDisc::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooOsipDisc

//______________________________________________________________________________
void rarStep::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarStep.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarStep(void *p) {
      return  p ? new(p) ::rarStep : new ::rarStep;
   }
   static void *newArray_rarStep(Long_t nElements, void *p) {
      return p ? new(p) ::rarStep[nElements] : new ::rarStep[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarStep(void *p) {
      delete ((::rarStep*)p);
   }
   static void deleteArray_rarStep(void *p) {
      delete [] ((::rarStep*)p);
   }
   static void destruct_rarStep(void *p) {
      typedef ::rarStep current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarStep(TBuffer &buf, void *obj) {
      ((::rarStep*)obj)->::rarStep::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarStep

//______________________________________________________________________________
void rarProd::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarProd.

   rarCompBase::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarProd(void *p) {
      return  p ? new(p) ::rarProd : new ::rarProd;
   }
   static void *newArray_rarProd(Long_t nElements, void *p) {
      return p ? new(p) ::rarProd[nElements] : new ::rarProd[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarProd(void *p) {
      delete ((::rarProd*)p);
   }
   static void deleteArray_rarProd(void *p) {
      delete [] ((::rarProd*)p);
   }
   static void destruct_rarProd(void *p) {
      typedef ::rarProd current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarProd(TBuffer &buf, void *obj) {
      ((::rarProd*)obj)->::rarProd::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarProd

//______________________________________________________________________________
void rarNLL::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarNLL.

   TNamed::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarNLL(void *p) {
      return  p ? new(p) ::rarNLL : new ::rarNLL;
   }
   static void *newArray_rarNLL(Long_t nElements, void *p) {
      return p ? new(p) ::rarNLL[nElements] : new ::rarNLL[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarNLL(void *p) {
      delete ((::rarNLL*)p);
   }
   static void deleteArray_rarNLL(void *p) {
      delete [] ((::rarNLL*)p);
   }
   static void destruct_rarNLL(void *p) {
      typedef ::rarNLL current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarNLL(TBuffer &buf, void *obj) {
      ((::rarNLL*)obj)->::rarNLL::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarNLL

//______________________________________________________________________________
void rarUniform::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarUniform.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarUniform(void *p) {
      return  p ? new(p) ::rarUniform : new ::rarUniform;
   }
   static void *newArray_rarUniform(Long_t nElements, void *p) {
      return p ? new(p) ::rarUniform[nElements] : new ::rarUniform[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarUniform(void *p) {
      delete ((::rarUniform*)p);
   }
   static void deleteArray_rarUniform(void *p) {
      delete [] ((::rarUniform*)p);
   }
   static void destruct_rarUniform(void *p) {
      typedef ::rarUniform current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarUniform(TBuffer &buf, void *obj) {
      ((::rarUniform*)obj)->::rarUniform::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarUniform

//______________________________________________________________________________
void rarOsipDisc::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarOsipDisc.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarOsipDisc(void *p) {
      return  p ? new(p) ::rarOsipDisc : new ::rarOsipDisc;
   }
   static void *newArray_rarOsipDisc(Long_t nElements, void *p) {
      return p ? new(p) ::rarOsipDisc[nElements] : new ::rarOsipDisc[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarOsipDisc(void *p) {
      delete ((::rarOsipDisc*)p);
   }
   static void deleteArray_rarOsipDisc(void *p) {
      delete [] ((::rarOsipDisc*)p);
   }
   static void destruct_rarOsipDisc(void *p) {
      typedef ::rarOsipDisc current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarOsipDisc(TBuffer &buf, void *obj) {
      ((::rarOsipDisc*)obj)->::rarOsipDisc::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarOsipDisc

//______________________________________________________________________________
void RooBallack::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooBallack.

   RooAbsPdf::Streamer(R__b);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RooBallack(void *p) {
      delete ((::RooBallack*)p);
   }
   static void deleteArray_RooBallack(void *p) {
      delete [] ((::RooBallack*)p);
   }
   static void destruct_RooBallack(void *p) {
      typedef ::RooBallack current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooBallack(TBuffer &buf, void *obj) {
      ((::RooBallack*)obj)->::RooBallack::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooBallack

//______________________________________________________________________________
void rarLass::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarLass.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarLass(void *p) {
      return  p ? new(p) ::rarLass : new ::rarLass;
   }
   static void *newArray_rarLass(Long_t nElements, void *p) {
      return p ? new(p) ::rarLass[nElements] : new ::rarLass[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarLass(void *p) {
      delete ((::rarLass*)p);
   }
   static void deleteArray_rarLass(void *p) {
      delete [] ((::rarLass*)p);
   }
   static void destruct_rarLass(void *p) {
      typedef ::rarLass current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarLass(TBuffer &buf, void *obj) {
      ((::rarLass*)obj)->::rarLass::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarLass

//______________________________________________________________________________
void rarAdd::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarAdd.

   rarCompBase::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarAdd(void *p) {
      return  p ? new(p) ::rarAdd : new ::rarAdd;
   }
   static void *newArray_rarAdd(Long_t nElements, void *p) {
      return p ? new(p) ::rarAdd[nElements] : new ::rarAdd[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarAdd(void *p) {
      delete ((::rarAdd*)p);
   }
   static void deleteArray_rarAdd(void *p) {
      delete [] ((::rarAdd*)p);
   }
   static void destruct_rarAdd(void *p) {
      typedef ::rarAdd current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarAdd(TBuffer &buf, void *obj) {
      ((::rarAdd*)obj)->::rarAdd::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarAdd

//______________________________________________________________________________
void rarTriGauss::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarTriGauss.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarTriGauss(void *p) {
      return  p ? new(p) ::rarTriGauss : new ::rarTriGauss;
   }
   static void *newArray_rarTriGauss(Long_t nElements, void *p) {
      return p ? new(p) ::rarTriGauss[nElements] : new ::rarTriGauss[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarTriGauss(void *p) {
      delete ((::rarTriGauss*)p);
   }
   static void deleteArray_rarTriGauss(void *p) {
      delete [] ((::rarTriGauss*)p);
   }
   static void destruct_rarTriGauss(void *p) {
      typedef ::rarTriGauss current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarTriGauss(TBuffer &buf, void *obj) {
      ((::rarTriGauss*)obj)->::rarTriGauss::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarTriGauss

//______________________________________________________________________________
void rarNovosibirsk::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarNovosibirsk.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarNovosibirsk(void *p) {
      return  p ? new(p) ::rarNovosibirsk : new ::rarNovosibirsk;
   }
   static void *newArray_rarNovosibirsk(Long_t nElements, void *p) {
      return p ? new(p) ::rarNovosibirsk[nElements] : new ::rarNovosibirsk[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarNovosibirsk(void *p) {
      delete ((::rarNovosibirsk*)p);
   }
   static void deleteArray_rarNovosibirsk(void *p) {
      delete [] ((::rarNovosibirsk*)p);
   }
   static void destruct_rarNovosibirsk(void *p) {
      typedef ::rarNovosibirsk current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarNovosibirsk(TBuffer &buf, void *obj) {
      ((::rarNovosibirsk*)obj)->::rarNovosibirsk::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarNovosibirsk

//______________________________________________________________________________
void rarHistPdf::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarHistPdf.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarHistPdf(void *p) {
      return  p ? new(p) ::rarHistPdf : new ::rarHistPdf;
   }
   static void *newArray_rarHistPdf(Long_t nElements, void *p) {
      return p ? new(p) ::rarHistPdf[nElements] : new ::rarHistPdf[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarHistPdf(void *p) {
      delete ((::rarHistPdf*)p);
   }
   static void deleteArray_rarHistPdf(void *p) {
      delete [] ((::rarHistPdf*)p);
   }
   static void destruct_rarHistPdf(void *p) {
      typedef ::rarHistPdf current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarHistPdf(TBuffer &buf, void *obj) {
      ((::rarHistPdf*)obj)->::rarHistPdf::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarHistPdf

//______________________________________________________________________________
void rarSPlot::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarSPlot.

   TNamed::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarSPlot(void *p) {
      return  p ? new(p) ::rarSPlot : new ::rarSPlot;
   }
   static void *newArray_rarSPlot(Long_t nElements, void *p) {
      return p ? new(p) ::rarSPlot[nElements] : new ::rarSPlot[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarSPlot(void *p) {
      delete ((::rarSPlot*)p);
   }
   static void deleteArray_rarSPlot(void *p) {
      delete [] ((::rarSPlot*)p);
   }
   static void destruct_rarSPlot(void *p) {
      typedef ::rarSPlot current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarSPlot(TBuffer &buf, void *obj) {
      ((::rarSPlot*)obj)->::rarSPlot::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarSPlot

//______________________________________________________________________________
void RooThreshold::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooThreshold.

   RooAbsPdf::Streamer(R__b);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RooThreshold(void *p) {
      delete ((::RooThreshold*)p);
   }
   static void deleteArray_RooThreshold(void *p) {
      delete [] ((::RooThreshold*)p);
   }
   static void destruct_RooThreshold(void *p) {
      typedef ::RooThreshold current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooThreshold(TBuffer &buf, void *obj) {
      ((::RooThreshold*)obj)->::RooThreshold::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooThreshold

//______________________________________________________________________________
void rarGaussian::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarGaussian.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarGaussian(void *p) {
      return  p ? new(p) ::rarGaussian : new ::rarGaussian;
   }
   static void *newArray_rarGaussian(Long_t nElements, void *p) {
      return p ? new(p) ::rarGaussian[nElements] : new ::rarGaussian[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarGaussian(void *p) {
      delete ((::rarGaussian*)p);
   }
   static void deleteArray_rarGaussian(void *p) {
      delete [] ((::rarGaussian*)p);
   }
   static void destruct_rarGaussian(void *p) {
      typedef ::rarGaussian current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarGaussian(TBuffer &buf, void *obj) {
      ((::rarGaussian*)obj)->::rarGaussian::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarGaussian

//______________________________________________________________________________
void RooCruijff::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooCruijff.

   RooAbsPdf::Streamer(R__b);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RooCruijff(void *p) {
      delete ((::RooCruijff*)p);
   }
   static void deleteArray_RooCruijff(void *p) {
      delete [] ((::RooCruijff*)p);
   }
   static void destruct_RooCruijff(void *p) {
      typedef ::RooCruijff current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooCruijff(TBuffer &buf, void *obj) {
      ((::RooCruijff*)obj)->::RooCruijff::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooCruijff

//______________________________________________________________________________
void rarBinned::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarBinned.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarBinned(void *p) {
      return  p ? new(p) ::rarBinned : new ::rarBinned;
   }
   static void *newArray_rarBinned(Long_t nElements, void *p) {
      return p ? new(p) ::rarBinned[nElements] : new ::rarBinned[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarBinned(void *p) {
      delete ((::rarBinned*)p);
   }
   static void deleteArray_rarBinned(void *p) {
      delete [] ((::rarBinned*)p);
   }
   static void destruct_rarBinned(void *p) {
      typedef ::rarBinned current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarBinned(TBuffer &buf, void *obj) {
      ((::rarBinned*)obj)->::rarBinned::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarBinned

//______________________________________________________________________________
void rarBallack::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarBallack.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarBallack(void *p) {
      return  p ? new(p) ::rarBallack : new ::rarBallack;
   }
   static void *newArray_rarBallack(Long_t nElements, void *p) {
      return p ? new(p) ::rarBallack[nElements] : new ::rarBallack[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarBallack(void *p) {
      delete ((::rarBallack*)p);
   }
   static void deleteArray_rarBallack(void *p) {
      delete [] ((::rarBallack*)p);
   }
   static void destruct_rarBallack(void *p) {
      typedef ::rarBallack current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarBallack(TBuffer &buf, void *obj) {
      ((::rarBallack*)obj)->::rarBallack::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarBallack

//______________________________________________________________________________
void rarRelBreitWigner::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarRelBreitWigner.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarRelBreitWigner(void *p) {
      return  p ? new(p) ::rarRelBreitWigner : new ::rarRelBreitWigner;
   }
   static void *newArray_rarRelBreitWigner(Long_t nElements, void *p) {
      return p ? new(p) ::rarRelBreitWigner[nElements] : new ::rarRelBreitWigner[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarRelBreitWigner(void *p) {
      delete ((::rarRelBreitWigner*)p);
   }
   static void deleteArray_rarRelBreitWigner(void *p) {
      delete [] ((::rarRelBreitWigner*)p);
   }
   static void destruct_rarRelBreitWigner(void *p) {
      typedef ::rarRelBreitWigner current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarRelBreitWigner(TBuffer &buf, void *obj) {
      ((::rarRelBreitWigner*)obj)->::rarRelBreitWigner::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarRelBreitWigner

//______________________________________________________________________________
void rarBifurGauss::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarBifurGauss.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarBifurGauss(void *p) {
      return  p ? new(p) ::rarBifurGauss : new ::rarBifurGauss;
   }
   static void *newArray_rarBifurGauss(Long_t nElements, void *p) {
      return p ? new(p) ::rarBifurGauss[nElements] : new ::rarBifurGauss[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarBifurGauss(void *p) {
      delete ((::rarBifurGauss*)p);
   }
   static void deleteArray_rarBifurGauss(void *p) {
      delete [] ((::rarBifurGauss*)p);
   }
   static void destruct_rarBifurGauss(void *p) {
      typedef ::rarBifurGauss current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarBifurGauss(TBuffer &buf, void *obj) {
      ((::rarBifurGauss*)obj)->::rarBifurGauss::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarBifurGauss

//______________________________________________________________________________
void rarUsrPdf::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarUsrPdf.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarUsrPdf(void *p) {
      return  p ? new(p) ::rarUsrPdf : new ::rarUsrPdf;
   }
   static void *newArray_rarUsrPdf(Long_t nElements, void *p) {
      return p ? new(p) ::rarUsrPdf[nElements] : new ::rarUsrPdf[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarUsrPdf(void *p) {
      delete ((::rarUsrPdf*)p);
   }
   static void deleteArray_rarUsrPdf(void *p) {
      delete [] ((::rarUsrPdf*)p);
   }
   static void destruct_rarUsrPdf(void *p) {
      typedef ::rarUsrPdf current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarUsrPdf(TBuffer &buf, void *obj) {
      ((::rarUsrPdf*)obj)->::rarUsrPdf::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarUsrPdf

//______________________________________________________________________________
void RooGounarisSakurai::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooGounarisSakurai.

   RooAbsPdf::Streamer(R__b);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RooGounarisSakurai(void *p) {
      delete ((::RooGounarisSakurai*)p);
   }
   static void deleteArray_RooGounarisSakurai(void *p) {
      delete [] ((::RooGounarisSakurai*)p);
   }
   static void destruct_RooGounarisSakurai(void *p) {
      typedef ::RooGounarisSakurai current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooGounarisSakurai(TBuffer &buf, void *obj) {
      ((::RooGounarisSakurai*)obj)->::RooGounarisSakurai::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooGounarisSakurai

//______________________________________________________________________________
void rarCBShape::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarCBShape.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarCBShape(void *p) {
      return  p ? new(p) ::rarCBShape : new ::rarCBShape;
   }
   static void *newArray_rarCBShape(Long_t nElements, void *p) {
      return p ? new(p) ::rarCBShape[nElements] : new ::rarCBShape[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarCBShape(void *p) {
      delete ((::rarCBShape*)p);
   }
   static void deleteArray_rarCBShape(void *p) {
      delete [] ((::rarCBShape*)p);
   }
   static void destruct_rarCBShape(void *p) {
      typedef ::rarCBShape current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarCBShape(TBuffer &buf, void *obj) {
      ((::rarCBShape*)obj)->::rarCBShape::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarCBShape

//______________________________________________________________________________
void rarGeneric::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarGeneric.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarGeneric(void *p) {
      return  p ? new(p) ::rarGeneric : new ::rarGeneric;
   }
   static void *newArray_rarGeneric(Long_t nElements, void *p) {
      return p ? new(p) ::rarGeneric[nElements] : new ::rarGeneric[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarGeneric(void *p) {
      delete ((::rarGeneric*)p);
   }
   static void deleteArray_rarGeneric(void *p) {
      delete [] ((::rarGeneric*)p);
   }
   static void destruct_rarGeneric(void *p) {
      typedef ::rarGeneric current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarGeneric(TBuffer &buf, void *obj) {
      ((::rarGeneric*)obj)->::rarGeneric::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarGeneric

//______________________________________________________________________________
void RooBinnedPdf::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooBinnedPdf.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      _x.Streamer(R__b);
      _coefList.Streamer(R__b);
      _limits.Streamer(R__b);
      R__b >> _nBins;
      R__b.CheckByteCount(R__s, R__c, RooBinnedPdf::IsA());
   } else {
      R__c = R__b.WriteVersion(RooBinnedPdf::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      _x.Streamer(R__b);
      _coefList.Streamer(R__b);
      _limits.Streamer(R__b);
      R__b << _nBins;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RooBinnedPdf(void *p) {
      delete ((::RooBinnedPdf*)p);
   }
   static void deleteArray_RooBinnedPdf(void *p) {
      delete [] ((::RooBinnedPdf*)p);
   }
   static void destruct_RooBinnedPdf(void *p) {
      typedef ::RooBinnedPdf current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooBinnedPdf(TBuffer &buf, void *obj) {
      ((::RooBinnedPdf*)obj)->::RooBinnedPdf::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooBinnedPdf

//______________________________________________________________________________
void rarCruijff::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarCruijff.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarCruijff(void *p) {
      return  p ? new(p) ::rarCruijff : new ::rarCruijff;
   }
   static void *newArray_rarCruijff(Long_t nElements, void *p) {
      return p ? new(p) ::rarCruijff[nElements] : new ::rarCruijff[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarCruijff(void *p) {
      delete ((::rarCruijff*)p);
   }
   static void deleteArray_rarCruijff(void *p) {
      delete [] ((::rarCruijff*)p);
   }
   static void destruct_rarCruijff(void *p) {
      typedef ::rarCruijff current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarCruijff(TBuffer &buf, void *obj) {
      ((::rarCruijff*)obj)->::rarCruijff::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarCruijff

//______________________________________________________________________________
void rarMinuit::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarMinuit.

   RooMinuit::Streamer(R__b);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_rarMinuit(void *p) {
      delete ((::rarMinuit*)p);
   }
   static void deleteArray_rarMinuit(void *p) {
      delete [] ((::rarMinuit*)p);
   }
   static void destruct_rarMinuit(void *p) {
      typedef ::rarMinuit current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarMinuit(TBuffer &buf, void *obj) {
      ((::rarMinuit*)obj)->::rarMinuit::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarMinuit

//______________________________________________________________________________
void rarSimPdf::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarSimPdf.

   rarCompBase::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarSimPdf(void *p) {
      return  p ? new(p) ::rarSimPdf : new ::rarSimPdf;
   }
   static void *newArray_rarSimPdf(Long_t nElements, void *p) {
      return p ? new(p) ::rarSimPdf[nElements] : new ::rarSimPdf[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarSimPdf(void *p) {
      delete ((::rarSimPdf*)p);
   }
   static void deleteArray_rarSimPdf(void *p) {
      delete [] ((::rarSimPdf*)p);
   }
   static void destruct_rarSimPdf(void *p) {
      typedef ::rarSimPdf current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarSimPdf(TBuffer &buf, void *obj) {
      ((::rarSimPdf*)obj)->::rarSimPdf::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarSimPdf

//______________________________________________________________________________
void rarFlatte::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarFlatte.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarFlatte(void *p) {
      return  p ? new(p) ::rarFlatte : new ::rarFlatte;
   }
   static void *newArray_rarFlatte(Long_t nElements, void *p) {
      return p ? new(p) ::rarFlatte[nElements] : new ::rarFlatte[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarFlatte(void *p) {
      delete ((::rarFlatte*)p);
   }
   static void deleteArray_rarFlatte(void *p) {
      delete [] ((::rarFlatte*)p);
   }
   static void destruct_rarFlatte(void *p) {
      typedef ::rarFlatte current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarFlatte(TBuffer &buf, void *obj) {
      ((::rarFlatte*)obj)->::rarFlatte::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarFlatte

//______________________________________________________________________________
void rarVoigtian::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarVoigtian.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarVoigtian(void *p) {
      return  p ? new(p) ::rarVoigtian : new ::rarVoigtian;
   }
   static void *newArray_rarVoigtian(Long_t nElements, void *p) {
      return p ? new(p) ::rarVoigtian[nElements] : new ::rarVoigtian[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarVoigtian(void *p) {
      delete ((::rarVoigtian*)p);
   }
   static void deleteArray_rarVoigtian(void *p) {
      delete [] ((::rarVoigtian*)p);
   }
   static void destruct_rarVoigtian(void *p) {
      typedef ::rarVoigtian current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarVoigtian(TBuffer &buf, void *obj) {
      ((::rarVoigtian*)obj)->::rarVoigtian::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarVoigtian

//______________________________________________________________________________
void rarThreshold::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarThreshold.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarThreshold(void *p) {
      return  p ? new(p) ::rarThreshold : new ::rarThreshold;
   }
   static void *newArray_rarThreshold(Long_t nElements, void *p) {
      return p ? new(p) ::rarThreshold[nElements] : new ::rarThreshold[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarThreshold(void *p) {
      delete ((::rarThreshold*)p);
   }
   static void deleteArray_rarThreshold(void *p) {
      delete [] ((::rarThreshold*)p);
   }
   static void destruct_rarThreshold(void *p) {
      typedef ::rarThreshold current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarThreshold(TBuffer &buf, void *obj) {
      ((::rarThreshold*)obj)->::rarThreshold::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarThreshold

//______________________________________________________________________________
void rarGaussModel::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarGaussModel.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarGaussModel(void *p) {
      return  p ? new(p) ::rarGaussModel : new ::rarGaussModel;
   }
   static void *newArray_rarGaussModel(Long_t nElements, void *p) {
      return p ? new(p) ::rarGaussModel[nElements] : new ::rarGaussModel[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarGaussModel(void *p) {
      delete ((::rarGaussModel*)p);
   }
   static void deleteArray_rarGaussModel(void *p) {
      delete [] ((::rarGaussModel*)p);
   }
   static void destruct_rarGaussModel(void *p) {
      typedef ::rarGaussModel current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarGaussModel(TBuffer &buf, void *obj) {
      ((::rarGaussModel*)obj)->::rarGaussModel::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarGaussModel

//______________________________________________________________________________
void rarExp::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarExp.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarExp(void *p) {
      return  p ? new(p) ::rarExp : new ::rarExp;
   }
   static void *newArray_rarExp(Long_t nElements, void *p) {
      return p ? new(p) ::rarExp[nElements] : new ::rarExp[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarExp(void *p) {
      delete ((::rarExp*)p);
   }
   static void deleteArray_rarExp(void *p) {
      delete [] ((::rarExp*)p);
   }
   static void destruct_rarExp(void *p) {
      typedef ::rarExp current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarExp(TBuffer &buf, void *obj) {
      ((::rarExp*)obj)->::rarExp::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarExp

//______________________________________________________________________________
void rarToyList::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarToyList.

   ::Error("rarToyList::Streamer", "version id <=0 in ClassDef, dummy Streamer() called"); if (R__b.IsReading()) { }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarToyList(void *p) {
      return  p ? new(p) ::rarToyList : new ::rarToyList;
   }
   static void *newArray_rarToyList(Long_t nElements, void *p) {
      return p ? new(p) ::rarToyList[nElements] : new ::rarToyList[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarToyList(void *p) {
      delete ((::rarToyList*)p);
   }
   static void deleteArray_rarToyList(void *p) {
      delete [] ((::rarToyList*)p);
   }
   static void destruct_rarToyList(void *p) {
      typedef ::rarToyList current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarToyList(TBuffer &buf, void *obj) {
      ((::rarToyList*)obj)->::rarToyList::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarToyList

//______________________________________________________________________________
void rarArgusBG::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarArgusBG.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarArgusBG(void *p) {
      return  p ? new(p) ::rarArgusBG : new ::rarArgusBG;
   }
   static void *newArray_rarArgusBG(Long_t nElements, void *p) {
      return p ? new(p) ::rarArgusBG[nElements] : new ::rarArgusBG[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarArgusBG(void *p) {
      delete ((::rarArgusBG*)p);
   }
   static void deleteArray_rarArgusBG(void *p) {
      delete [] ((::rarArgusBG*)p);
   }
   static void destruct_rarArgusBG(void *p) {
      typedef ::rarArgusBG current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarArgusBG(TBuffer &buf, void *obj) {
      ((::rarArgusBG*)obj)->::rarArgusBG::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarArgusBG

//______________________________________________________________________________
void RooRelBreitWigner::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooRelBreitWigner.

   RooAbsPdf::Streamer(R__b);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RooRelBreitWigner(void *p) {
      delete ((::RooRelBreitWigner*)p);
   }
   static void deleteArray_RooRelBreitWigner(void *p) {
      delete [] ((::RooRelBreitWigner*)p);
   }
   static void destruct_RooRelBreitWigner(void *p) {
      typedef ::RooRelBreitWigner current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooRelBreitWigner(TBuffer &buf, void *obj) {
      ((::RooRelBreitWigner*)obj)->::RooRelBreitWigner::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooRelBreitWigner

//______________________________________________________________________________
void rarGounarisSakurai::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarGounarisSakurai.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarGounarisSakurai(void *p) {
      return  p ? new(p) ::rarGounarisSakurai : new ::rarGounarisSakurai;
   }
   static void *newArray_rarGounarisSakurai(Long_t nElements, void *p) {
      return p ? new(p) ::rarGounarisSakurai[nElements] : new ::rarGounarisSakurai[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarGounarisSakurai(void *p) {
      delete ((::rarGounarisSakurai*)p);
   }
   static void deleteArray_rarGounarisSakurai(void *p) {
      delete [] ((::rarGounarisSakurai*)p);
   }
   static void destruct_rarGounarisSakurai(void *p) {
      typedef ::rarGounarisSakurai current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarGounarisSakurai(TBuffer &buf, void *obj) {
      ((::rarGounarisSakurai*)obj)->::rarGounarisSakurai::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarGounarisSakurai

//______________________________________________________________________________
void rarMLFitter::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarMLFitter.

   rarCompBase::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarMLFitter(void *p) {
      return  p ? new(p) ::rarMLFitter : new ::rarMLFitter;
   }
   static void *newArray_rarMLFitter(Long_t nElements, void *p) {
      return p ? new(p) ::rarMLFitter[nElements] : new ::rarMLFitter[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarMLFitter(void *p) {
      delete ((::rarMLFitter*)p);
   }
   static void deleteArray_rarMLFitter(void *p) {
      delete [] ((::rarMLFitter*)p);
   }
   static void destruct_rarMLFitter(void *p) {
      typedef ::rarMLFitter current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarMLFitter(TBuffer &buf, void *obj) {
      ((::rarMLFitter*)obj)->::rarMLFitter::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarMLFitter

//______________________________________________________________________________
void rarKeys::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarKeys.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarKeys(void *p) {
      return  p ? new(p) ::rarKeys : new ::rarKeys;
   }
   static void *newArray_rarKeys(Long_t nElements, void *p) {
      return p ? new(p) ::rarKeys[nElements] : new ::rarKeys[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarKeys(void *p) {
      delete ((::rarKeys*)p);
   }
   static void deleteArray_rarKeys(void *p) {
      delete [] ((::rarKeys*)p);
   }
   static void destruct_rarKeys(void *p) {
      typedef ::rarKeys current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarKeys(TBuffer &buf, void *obj) {
      ((::rarKeys*)obj)->::rarKeys::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarKeys

//______________________________________________________________________________
void rarMLPdf::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarMLPdf.

   rarAdd::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarMLPdf(void *p) {
      return  p ? new(p) ::rarMLPdf : new ::rarMLPdf;
   }
   static void *newArray_rarMLPdf(Long_t nElements, void *p) {
      return p ? new(p) ::rarMLPdf[nElements] : new ::rarMLPdf[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarMLPdf(void *p) {
      delete ((::rarMLPdf*)p);
   }
   static void deleteArray_rarMLPdf(void *p) {
      delete [] ((::rarMLPdf*)p);
   }
   static void destruct_rarMLPdf(void *p) {
      typedef ::rarMLPdf current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarMLPdf(TBuffer &buf, void *obj) {
      ((::rarMLPdf*)obj)->::rarMLPdf::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarMLPdf

//______________________________________________________________________________
void RooFlatte::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooFlatte.

   RooAbsPdf::Streamer(R__b);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RooFlatte(void *p) {
      delete ((::RooFlatte*)p);
   }
   static void deleteArray_RooFlatte(void *p) {
      delete [] ((::RooFlatte*)p);
   }
   static void destruct_RooFlatte(void *p) {
      typedef ::RooFlatte current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooFlatte(TBuffer &buf, void *obj) {
      ((::RooFlatte*)obj)->::RooFlatte::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooFlatte

//______________________________________________________________________________
void rarDecay::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarDecay.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarDecay(void *p) {
      return  p ? new(p) ::rarDecay : new ::rarDecay;
   }
   static void *newArray_rarDecay(Long_t nElements, void *p) {
      return p ? new(p) ::rarDecay[nElements] : new ::rarDecay[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarDecay(void *p) {
      delete ((::rarDecay*)p);
   }
   static void deleteArray_rarDecay(void *p) {
      delete [] ((::rarDecay*)p);
   }
   static void destruct_rarDecay(void *p) {
      typedef ::rarDecay current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarDecay(TBuffer &buf, void *obj) {
      ((::rarDecay*)obj)->::rarDecay::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarDecay

//______________________________________________________________________________
void rarPoly::Streamer(TBuffer &R__b)
{
   // Stream an object of class rarPoly.

   rarBasePdf::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_rarPoly(void *p) {
      return  p ? new(p) ::rarPoly : new ::rarPoly;
   }
   static void *newArray_rarPoly(Long_t nElements, void *p) {
      return p ? new(p) ::rarPoly[nElements] : new ::rarPoly[nElements];
   }
   // Wrapper around operator delete
   static void delete_rarPoly(void *p) {
      delete ((::rarPoly*)p);
   }
   static void deleteArray_rarPoly(void *p) {
      delete [] ((::rarPoly*)p);
   }
   static void destruct_rarPoly(void *p) {
      typedef ::rarPoly current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_rarPoly(TBuffer &buf, void *obj) {
      ((::rarPoly*)obj)->::rarPoly::Streamer(buf);
   }
} // end of namespace ROOT for class ::rarPoly

namespace {
  void TriggerDictionaryInitialization_RooRarFitCint_Impl() {
    static const char* headers[] = {
"rarStrParser.hh",
"RooLass.hh",
"rarTwoGauss.hh",
"rarMultPdf.hh",
"RooOsipDisc.hh",
"rarStep.hh",
"rarProd.hh",
"rarNLL.hh",
"rarUniform.hh",
"rarConfig.hh",
"rarOsipDisc.hh",
"RooBallack.hh",
"rarLass.hh",
"rarAdd.hh",
"rarTriGauss.hh",
"rarNovosibirsk.hh",
"rarHistPdf.hh",
"rarSPlot.hh",
"RooThreshold.hh",
"rarGaussian.hh",
"RooCruijff.hh",
"rarBasePdf.hh",
"rarDatasets.hh",
"rarBinned.hh",
"rarBallack.hh",
"rarRelBreitWigner.hh",
"rarBifurGauss.hh",
"rarUsrPdf.hh",
"rarDatasetDef.hh",
"RooGounarisSakurai.hh",
"rarCompBase.hh",
"rarCBShape.hh",
"rarGeneric.hh",
"RooBinnedPdf.hh",
"rarCruijff.hh",
"rarMinuit.hh",
"rarSimPdf.hh",
"rarFlatte.hh",
"rarVoigtian.hh",
"rarThreshold.hh",
"rarGaussModel.hh",
"rarExp.hh",
"rarToyList.hh",
"rarArgusBG.hh",
"RooRelBreitWigner.hh",
"rarGounarisSakurai.hh",
"rarMLFitter.hh",
"rarKeys.hh",
"rarMLPdf.hh",
"RooFlatte.hh",
"rarDecay.hh",
"rarPoly.hh",
0
    };
    static const char* includePaths[] = {
"..",
"/pbs/home/a/azakaria/rich_matrices/RooRarFit/tmp/RooRarFit/",
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.14.04-4d676/x86_64-centos7-gcc62-opt/include",
"/pbs/home/a/azakaria/rich_matrices/RooRarFit/tmp/RooRarFit/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "RooRarFitCint dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit String Parser class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarStrParser.hh")))  rarStrParser;
class __attribute__((annotate("$clingAutoload$RooLass.hh")))  RooLass;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit config class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarConfig.hh")))  __attribute__((annotate("$clingAutoload$rarTwoGauss.hh")))  rarConfig;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit dataset definition class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarDatasetDef.hh")))  __attribute__((annotate("$clingAutoload$rarTwoGauss.hh")))  rarDatasetDef;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit dataset class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarDatasets.hh")))  __attribute__((annotate("$clingAutoload$rarTwoGauss.hh")))  rarDatasets;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit base Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarBasePdf.hh")))  __attribute__((annotate("$clingAutoload$rarTwoGauss.hh")))  rarBasePdf;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit TwoGauss Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarTwoGauss.hh")))  rarTwoGauss;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Component base class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarCompBase.hh")))  __attribute__((annotate("$clingAutoload$rarMultPdf.hh")))  rarCompBase;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Multiple class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarMultPdf.hh")))  rarMultPdf;
class __attribute__((annotate("$clingAutoload$RooOsipDisc.hh")))  RooOsipDisc;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit ParametricStep Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarStep.hh")))  rarStep;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit ProdPdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarProd.hh")))  rarProd;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit NLL class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarNLL.hh")))  rarNLL;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Uniform Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarUniform.hh")))  rarUniform;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Osipenkov pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarOsipDisc.hh")))  rarOsipDisc;
class __attribute__((annotate("$clingAutoload$RooBallack.hh")))  RooBallack;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit User-defined Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarLass.hh")))  rarLass;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit AddPdf/AddModel class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarAdd.hh")))  rarAdd;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit TriGauss related Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarTriGauss.hh")))  rarTriGauss;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Novosibirsk Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarNovosibirsk.hh")))  rarNovosibirsk;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit HistPdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarHistPdf.hh")))  rarHistPdf;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit sPlot class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarSPlot.hh")))  rarSPlot;
class __attribute__((annotate(R"ATTRDUMP(Threshold PDF)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RooThreshold.hh")))  RooThreshold;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Gaussian/BreitWigner Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarGaussian.hh")))  rarGaussian;
class __attribute__((annotate("$clingAutoload$RooCruijff.hh")))  RooCruijff;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Binned Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarBinned.hh")))  rarBinned;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Ballack Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarBallack.hh")))  rarBallack;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit RelBreitWigner PDF class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarRelBreitWigner.hh")))  rarRelBreitWigner;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit BifurGauss Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarBifurGauss.hh")))  rarBifurGauss;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit User-defined Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarUsrPdf.hh")))  rarUsrPdf;
class __attribute__((annotate(R"ATTRDUMP(Gounaris Sakurai PDF)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RooGounarisSakurai.hh")))  RooGounarisSakurai;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit CBShape Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarCBShape.hh")))  rarCBShape;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Generic Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarGeneric.hh")))  rarGeneric;
class __attribute__((annotate(R"ATTRDUMP(Parametric Step Function Pdf)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RooBinnedPdf.hh")))  RooBinnedPdf;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Cruijff Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarCruijff.hh")))  rarCruijff;
class __attribute__((annotate(R"ATTRDUMP(RooMinuit derivative with contour RooPlot)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarMinuit.hh")))  rarMinuit;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Simultaneous Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarSimPdf.hh")))  rarSimPdf;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Flatte PDF class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarFlatte.hh")))  rarFlatte;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Voigtian Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarVoigtian.hh")))  rarVoigtian;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Threshold PDF class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarThreshold.hh")))  rarThreshold;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Gaussian Resolution Model class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarGaussModel.hh")))  rarGaussModel;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Exponential Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarExp.hh")))  rarExp;
class __attribute__((annotate("$clingAutoload$rarToyList.hh")))  rarToyList;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit ArgusBG Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarArgusBG.hh")))  rarArgusBG;
class __attribute__((annotate(R"ATTRDUMP(Relativistic Breit Wigner PDF)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RooRelBreitWigner.hh")))  RooRelBreitWigner;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit GounarisSakurai PDF class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarGounarisSakurai.hh")))  rarGounarisSakurai;
class __attribute__((annotate(R"ATTRDUMP(Final RooRarFit ML fitter class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarMLFitter.hh")))  rarMLFitter;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit 1D/2D Keys Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarKeys.hh")))  rarKeys;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit ML model class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarMLPdf.hh")))  rarMLPdf;
class __attribute__((annotate(R"ATTRDUMP(Flatte PDF using BES parameterisation)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RooFlatte.hh")))  RooFlatte;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit BDecay/Decay model class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarDecay.hh")))  rarDecay;
class __attribute__((annotate(R"ATTRDUMP(RooRarFit Polynomial/Chebychev Pdf class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$rarPoly.hh")))  rarPoly;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "RooRarFitCint dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "rarStrParser.hh"
#include "RooLass.hh"
#include "rarTwoGauss.hh"
#include "rarMultPdf.hh"
#include "RooOsipDisc.hh"
#include "rarStep.hh"
#include "rarProd.hh"
#include "rarNLL.hh"
#include "rarUniform.hh"
#include "rarConfig.hh"
#include "rarOsipDisc.hh"
#include "RooBallack.hh"
#include "rarLass.hh"
#include "rarAdd.hh"
#include "rarTriGauss.hh"
#include "rarNovosibirsk.hh"
#include "rarHistPdf.hh"
#include "rarSPlot.hh"
#include "RooThreshold.hh"
#include "rarGaussian.hh"
#include "RooCruijff.hh"
#include "rarBasePdf.hh"
#include "rarDatasets.hh"
#include "rarBinned.hh"
#include "rarBallack.hh"
#include "rarRelBreitWigner.hh"
#include "rarBifurGauss.hh"
#include "rarUsrPdf.hh"
#include "rarDatasetDef.hh"
#include "RooGounarisSakurai.hh"
#include "rarCompBase.hh"
#include "rarCBShape.hh"
#include "rarGeneric.hh"
#include "RooBinnedPdf.hh"
#include "rarCruijff.hh"
#include "rarMinuit.hh"
#include "rarSimPdf.hh"
#include "rarFlatte.hh"
#include "rarVoigtian.hh"
#include "rarThreshold.hh"
#include "rarGaussModel.hh"
#include "rarExp.hh"
#include "rarToyList.hh"
#include "rarArgusBG.hh"
#include "RooRelBreitWigner.hh"
#include "rarGounarisSakurai.hh"
#include "rarMLFitter.hh"
#include "rarKeys.hh"
#include "rarMLPdf.hh"
#include "RooFlatte.hh"
#include "rarDecay.hh"
#include "rarPoly.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"RooBallack", payloadCode, "@",
"RooBinnedPdf", payloadCode, "@",
"RooCruijff", payloadCode, "@",
"RooFlatte", payloadCode, "@",
"RooGounarisSakurai", payloadCode, "@",
"RooLass", payloadCode, "@",
"RooOsipDisc", payloadCode, "@",
"RooRelBreitWigner", payloadCode, "@",
"RooThreshold", payloadCode, "@",
"rarAdd", payloadCode, "@",
"rarArgusBG", payloadCode, "@",
"rarBallack", payloadCode, "@",
"rarBasePdf", payloadCode, "@",
"rarBifurGauss", payloadCode, "@",
"rarBinned", payloadCode, "@",
"rarCBShape", payloadCode, "@",
"rarCompBase", payloadCode, "@",
"rarConfig", payloadCode, "@",
"rarCruijff", payloadCode, "@",
"rarDatasetDef", payloadCode, "@",
"rarDatasets", payloadCode, "@",
"rarDecay", payloadCode, "@",
"rarExp", payloadCode, "@",
"rarFlatte", payloadCode, "@",
"rarGaussModel", payloadCode, "@",
"rarGaussian", payloadCode, "@",
"rarGeneric", payloadCode, "@",
"rarGounarisSakurai", payloadCode, "@",
"rarHistPdf", payloadCode, "@",
"rarKeys", payloadCode, "@",
"rarLass", payloadCode, "@",
"rarMLFitter", payloadCode, "@",
"rarMLPdf", payloadCode, "@",
"rarMinuit", payloadCode, "@",
"rarMultPdf", payloadCode, "@",
"rarNLL", payloadCode, "@",
"rarNovosibirsk", payloadCode, "@",
"rarOsipDisc", payloadCode, "@",
"rarPoly", payloadCode, "@",
"rarProd", payloadCode, "@",
"rarRelBreitWigner", payloadCode, "@",
"rarSPlot", payloadCode, "@",
"rarSimPdf", payloadCode, "@",
"rarStep", payloadCode, "@",
"rarStrParser", payloadCode, "@",
"rarThreshold", payloadCode, "@",
"rarToyList", payloadCode, "@",
"rarTriGauss", payloadCode, "@",
"rarTwoGauss", payloadCode, "@",
"rarUniform", payloadCode, "@",
"rarUsrPdf", payloadCode, "@",
"rarVoigtian", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("RooRarFitCint",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_RooRarFitCint_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_RooRarFitCint_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_RooRarFitCint() {
  TriggerDictionaryInitialization_RooRarFitCint_Impl();
}
