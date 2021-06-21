// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIsrcdIHitsdict
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
#include "Hits.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_Hits(void *p = 0);
   static void *newArray_Hits(Long_t size, void *p);
   static void delete_Hits(void *p);
   static void deleteArray_Hits(void *p);
   static void destruct_Hits(void *p);
   static void streamer_Hits(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Hits*)
   {
      ::Hits *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Hits >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Hits", ::Hits::Class_Version(), "Hits.h", 38,
                  typeid(::Hits), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Hits::Dictionary, isa_proxy, 16,
                  sizeof(::Hits) );
      instance.SetNew(&new_Hits);
      instance.SetNewArray(&newArray_Hits);
      instance.SetDelete(&delete_Hits);
      instance.SetDeleteArray(&deleteArray_Hits);
      instance.SetDestructor(&destruct_Hits);
      instance.SetStreamerFunc(&streamer_Hits);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Hits*)
   {
      return GenerateInitInstanceLocal((::Hits*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Hits*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Hits::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Hits::Class_Name()
{
   return "Hits";
}

//______________________________________________________________________________
const char *Hits::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Hits*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Hits::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Hits*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Hits::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Hits*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Hits::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Hits*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Hits::Streamer(TBuffer &R__b)
{
   // Stream an object of class Hits.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> ENum;
      R__b >> NEle;
      R__b >> NHits;
      R__b >> NPrimHits;
      R__b >> NClus;
      R__b >> NTrips;
      R__b >> NFitUp;
      R__b >> NFitDown;
      R__b >> NFinders;
      R__b >> NParticles;
      R__b >> NRecTracks;
      R__b >> NShowerHits;
      fLastHit.Streamer(R__b);
      fHits->Streamer(R__b);
      R__b >> fIsValid;
      R__b.CheckByteCount(R__s, R__c, Hits::IsA());
   } else {
      R__c = R__b.WriteVersion(Hits::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << ENum;
      R__b << NEle;
      R__b << NHits;
      R__b << NPrimHits;
      R__b << NClus;
      R__b << NTrips;
      R__b << NFitUp;
      R__b << NFitDown;
      R__b << NFinders;
      R__b << NParticles;
      R__b << NRecTracks;
      R__b << NShowerHits;
      fLastHit.Streamer(R__b);
      fHits->Streamer(R__b);
      R__b << fIsValid;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Hits(void *p) {
      return  p ? new(p) ::Hits : new ::Hits;
   }
   static void *newArray_Hits(Long_t nElements, void *p) {
      return p ? new(p) ::Hits[nElements] : new ::Hits[nElements];
   }
   // Wrapper around operator delete
   static void delete_Hits(void *p) {
      delete ((::Hits*)p);
   }
   static void deleteArray_Hits(void *p) {
      delete [] ((::Hits*)p);
   }
   static void destruct_Hits(void *p) {
      typedef ::Hits current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_Hits(TBuffer &buf, void *obj) {
      ((::Hits*)obj)->::Hits::Streamer(buf);
   }
} // end of namespace ROOT for class ::Hits

namespace {
  void TriggerDictionaryInitialization_Hitsdict_Impl() {
    static const char* headers[] = {
"Hits.h",
0
    };
    static const char* includePaths[] = {
"/products/root/root-6.20.00/include/",
"/home/suryanarayan/Documents/Gobinda/IICHEP/mical_20190829/include/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "Hitsdict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Event structure)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$Hits.h")))  Hits;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "Hitsdict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "Hits.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"Hits", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("Hitsdict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_Hitsdict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_Hitsdict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_Hitsdict() {
  TriggerDictionaryInitialization_Hitsdict_Impl();
}
