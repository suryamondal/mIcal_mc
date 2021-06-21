// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIsrcdIHitPosdict
#define R__NO_DEPRECATION

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

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "HitPos.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_HitPos(void *p = 0);
   static void *newArray_HitPos(Long_t size, void *p);
   static void delete_HitPos(void *p);
   static void deleteArray_HitPos(void *p);
   static void destruct_HitPos(void *p);
   static void streamer_HitPos(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HitPos*)
   {
      ::HitPos *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HitPos >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HitPos", ::HitPos::Class_Version(), "HitPos.h", 34,
                  typeid(::HitPos), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HitPos::Dictionary, isa_proxy, 16,
                  sizeof(::HitPos) );
      instance.SetNew(&new_HitPos);
      instance.SetNewArray(&newArray_HitPos);
      instance.SetDelete(&delete_HitPos);
      instance.SetDeleteArray(&deleteArray_HitPos);
      instance.SetDestructor(&destruct_HitPos);
      instance.SetStreamerFunc(&streamer_HitPos);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HitPos*)
   {
      return GenerateInitInstanceLocal((::HitPos*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HitPos*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr HitPos::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HitPos::Class_Name()
{
   return "HitPos";
}

//______________________________________________________________________________
const char *HitPos::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HitPos*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HitPos::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HitPos*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HitPos::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HitPos*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HitPos::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HitPos*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void HitPos::Streamer(TBuffer &R__b)
{
   // Stream an object of class HitPos.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> TrackType;
      R__b >> HitNum;
      R__b >> PrimHitNum;
      R__b >> CluNum;
      R__b >> TriNum;
      R__b >> FitUpNum;
      R__b >> FitDownNum;
      R__b >> FindNum;
      R__b >> ShowerHitNum;
      R__b >> ParCode;
      R__b >> Fup;
      R__b >> XX;
      R__b >> YY;
      R__b >> ZZ;
      R__b >> pmag;
      R__b >> pt;
      R__b >> pp;
      R__b.CheckByteCount(R__s, R__c, HitPos::IsA());
   } else {
      R__c = R__b.WriteVersion(HitPos::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << TrackType;
      R__b << HitNum;
      R__b << PrimHitNum;
      R__b << CluNum;
      R__b << TriNum;
      R__b << FitUpNum;
      R__b << FitDownNum;
      R__b << FindNum;
      R__b << ShowerHitNum;
      R__b << ParCode;
      R__b << Fup;
      R__b << XX;
      R__b << YY;
      R__b << ZZ;
      R__b << pmag;
      R__b << pt;
      R__b << pp;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HitPos(void *p) {
      return  p ? new(p) ::HitPos : new ::HitPos;
   }
   static void *newArray_HitPos(Long_t nElements, void *p) {
      return p ? new(p) ::HitPos[nElements] : new ::HitPos[nElements];
   }
   // Wrapper around operator delete
   static void delete_HitPos(void *p) {
      delete ((::HitPos*)p);
   }
   static void deleteArray_HitPos(void *p) {
      delete [] ((::HitPos*)p);
   }
   static void destruct_HitPos(void *p) {
      typedef ::HitPos current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_HitPos(TBuffer &buf, void *obj) {
      ((::HitPos*)obj)->::HitPos::Streamer(buf);
   }
} // end of namespace ROOT for class ::HitPos

namespace {
  void TriggerDictionaryInitialization_HitPosdict_Impl() {
    static const char* headers[] = {
"HitPos.h",
0
    };
    static const char* includePaths[] = {
"/products/root/root-6.20.00/include/",
"/home/suryanarayan/Documents/Gobinda/IICHEP/mical_20190829/include/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "HitPosdict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(A track segment)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$HitPos.h")))  HitPos;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "HitPosdict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "HitPos.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"HitPos", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("HitPosdict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_HitPosdict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_HitPosdict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_HitPosdict() {
  TriggerDictionaryInitialization_HitPosdict_Impl();
}
