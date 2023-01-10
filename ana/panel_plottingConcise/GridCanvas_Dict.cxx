// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME GridCanvas_Dict

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
#include "GridCanvas.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_GridCanvas(void *p = 0);
   static void *newArray_GridCanvas(Long_t size, void *p);
   static void delete_GridCanvas(void *p);
   static void deleteArray_GridCanvas(void *p);
   static void destruct_GridCanvas(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::GridCanvas*)
   {
      ::GridCanvas *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::GridCanvas >(0);
      static ::ROOT::TGenericClassInfo 
         instance("GridCanvas", ::GridCanvas::Class_Version(), "GridCanvas.h", 13,
                  typeid(::GridCanvas), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::GridCanvas::Dictionary, isa_proxy, 4,
                  sizeof(::GridCanvas) );
      instance.SetNew(&new_GridCanvas);
      instance.SetNewArray(&newArray_GridCanvas);
      instance.SetDelete(&delete_GridCanvas);
      instance.SetDeleteArray(&deleteArray_GridCanvas);
      instance.SetDestructor(&destruct_GridCanvas);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::GridCanvas*)
   {
      return GenerateInitInstanceLocal((::GridCanvas*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::GridCanvas*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr GridCanvas::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *GridCanvas::Class_Name()
{
   return "GridCanvas";
}

//______________________________________________________________________________
const char *GridCanvas::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GridCanvas*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int GridCanvas::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GridCanvas*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GridCanvas::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GridCanvas*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GridCanvas::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GridCanvas*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void GridCanvas::Streamer(TBuffer &R__b)
{
   // Stream an object of class GridCanvas.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(GridCanvas::Class(),this);
   } else {
      R__b.WriteClassBuffer(GridCanvas::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_GridCanvas(void *p) {
      return  p ? new(p) ::GridCanvas : new ::GridCanvas;
   }
   static void *newArray_GridCanvas(Long_t nElements, void *p) {
      return p ? new(p) ::GridCanvas[nElements] : new ::GridCanvas[nElements];
   }
   // Wrapper around operator delete
   static void delete_GridCanvas(void *p) {
      delete ((::GridCanvas*)p);
   }
   static void deleteArray_GridCanvas(void *p) {
      delete [] ((::GridCanvas*)p);
   }
   static void destruct_GridCanvas(void *p) {
      typedef ::GridCanvas current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::GridCanvas

namespace {
  void TriggerDictionaryInitialization_GridCanvas_Dict_Impl() {
    static const char* headers[] = {
"GridCanvas.h",
0
    };
    static const char* includePaths[] = {
"/minerva/app/users/anezkak/MAT_GitHub/opt/lib/../include/",
"/cvmfs/minerva.opensciencegrid.org/minerva/hep_hpc_products/root/v6_10_04d/Linux64bit+3.10-2.17-e14-prof/include",
"/minerva/app/users/anezkak/MAT_GitHub/NSFNukeCCInclusive/ana/panel_plottingConcise/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "GridCanvas_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$GridCanvas.h")))  GridCanvas;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "GridCanvas_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "GridCanvas.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"GridCanvas", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("GridCanvas_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_GridCanvas_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_GridCanvas_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_GridCanvas_Dict() {
  TriggerDictionaryInitialization_GridCanvas_Dict_Impl();
}
