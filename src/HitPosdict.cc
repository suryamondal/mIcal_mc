//
// File generated by rootcint at Tue Dec 31 10:35:10 2019

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME dOdOdIsrcdIHitPosdict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "HitPosdict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// Direct notice to TROOT of the dictionary's loading.
namespace {
   static struct DictInit {
      DictInit() {
         ROOT::RegisterModule();
      }
   } __TheDictionaryInitializer;
}

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void HitPos_ShowMembers(void *obj, TMemberInspector &R__insp);
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
         instance("HitPos", ::HitPos::Class_Version(), "./HitPos.h", 34,
                  typeid(::HitPos), DefineBehavior(ptr, ptr),
                  &::HitPos::Dictionary, isa_proxy, 0,
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
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::HitPos*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *HitPos::fgIsA = 0;  // static to hold class pointer

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
void HitPos::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HitPos*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *HitPos::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HitPos*)0x0)->GetClass();
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

//______________________________________________________________________________
void HitPos::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class HitPos.
      TClass *R__cl = ::HitPos::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "TrackType", &TrackType);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "HitNum", &HitNum);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "PrimHitNum", &PrimHitNum);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "CluNum", &CluNum);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "TriNum", &TriNum);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "FitUpNum", &FitUpNum);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "FitDownNum", &FitDownNum);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "FindNum", &FindNum);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "ShowerHitNum", &ShowerHitNum);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "ParCode", &ParCode);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Fup", &Fup);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "XX", &XX);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "YY", &YY);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "ZZ", &ZZ);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "pmag", &pmag);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "pt", &pt);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "pp", &pp);
      TObject::ShowMembers(R__insp);
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

/********************************************************
* ../src/HitPosdict.cc
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableHitPosdict();

extern "C" void G__set_cpp_environmentHitPosdict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("HitPos.h");
  G__cpp_reset_tagtableHitPosdict();
}
#include <new>
extern "C" int G__cpp_dllrevHitPosdict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* HitPos */
static int G__HitPosdict_233_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   HitPos* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new HitPos[n];
     } else {
       p = new((void*) gvp) HitPos[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new HitPos;
     } else {
       p = new((void*) gvp) HitPos;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__HitPosdictLN_HitPos));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitPosdict_233_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   HitPos* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new HitPos((Float_t) G__double(libp->para[0]));
   } else {
     p = new((void*) gvp) HitPos((Float_t) G__double(libp->para[0]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__HitPosdictLN_HitPos));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitPosdict_233_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) HitPos::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitPosdict_233_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) HitPos::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitPosdict_233_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) HitPos::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitPosdict_233_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      HitPos::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitPosdict_233_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((HitPos*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitPosdict_233_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) HitPos::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitPosdict_233_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) HitPos::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitPosdict_233_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) HitPos::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HitPosdict_233_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) HitPos::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__HitPosdict_233_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   HitPos* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new HitPos(*(HitPos*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__HitPosdictLN_HitPos));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef HitPos G__THitPos;
static int G__HitPosdict_233_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (HitPos*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((HitPos*) (soff+(sizeof(HitPos)*i)))->~G__THitPos();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (HitPos*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((HitPos*) (soff))->~G__THitPos();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__HitPosdict_233_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   HitPos* dest = (HitPos*) G__getstructoffset();
   *dest = *(HitPos*) libp->para[0].ref;
   const HitPos& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* HitPos */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncHitPosdict {
 public:
  G__Sizep2memfuncHitPosdict(): p(&G__Sizep2memfuncHitPosdict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncHitPosdict::*p)();
};

size_t G__get_sizep2memfuncHitPosdict()
{
  G__Sizep2memfuncHitPosdict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceHitPosdict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__HitPosdictLN_HitPos))) {
     HitPos *G__Lderived;
     G__Lderived=(HitPos*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__HitPosdictLN_HitPos),G__get_linked_tagnum(&G__HitPosdictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableHitPosdict() {

   /* Setting up typedef entry */
   G__search_typename2("Float_t",102,-1,0,-1);
   G__setnewtype(-1,"Float 4 bytes (float)",0);
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__HitPosdictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__HitPosdictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__HitPosdictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__HitPosdictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__HitPosdictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__HitPosdictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__HitPosdictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__HitPosdictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__HitPosdictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__HitPosdictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<std::bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__HitPosdictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__HitPosdictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*>",117,G__get_linked_tagnum(&G__HitPosdictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long>",117,G__get_linked_tagnum(&G__HitPosdictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long,const TObject**>",117,G__get_linked_tagnum(&G__HitPosdictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TVectorT<Float_t>",117,G__get_linked_tagnum(&G__HitPosdictLN_TVectorTlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TVectorT<Double_t>",117,G__get_linked_tagnum(&G__HitPosdictLN_TVectorTlEdoublegR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* HitPos */
static void G__setup_memvarHitPos(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__HitPosdictLN_HitPos));
   { HitPos *p; p=(HitPos*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->TrackType)-(long)(p)),105,0,0,-1,-1,-1,1,"TrackType=",0,"Track Type: -1: hits, -101 : Primary hits -2: clulster, -3: triplet, -4: track,-44: track /fitter?");
   G__memvar_setup((void*)((long)(&p->HitNum)-(long)(p)),105,0,0,-1,-1,-1,1,"HitNum=",0,"Hit Number");
   G__memvar_setup((void*)((long)(&p->PrimHitNum)-(long)(p)),105,0,0,-1,-1,-1,1,"PrimHitNum=",0,"Primary Number ");
   G__memvar_setup((void*)((long)(&p->CluNum)-(long)(p)),105,0,0,-1,-1,-1,1,"CluNum=",0,"Cluster Number");
   G__memvar_setup((void*)((long)(&p->TriNum)-(long)(p)),105,0,0,-1,-1,-1,1,"TriNum=",0,"Triplet Number");
   G__memvar_setup((void*)((long)(&p->FitUpNum)-(long)(p)),105,0,0,-1,-1,-1,1,"FitUpNum=",0,"Fit  Number");
   G__memvar_setup((void*)((long)(&p->FitDownNum)-(long)(p)),105,0,0,-1,-1,-1,1,"FitDownNum=",0,"Fit  Number");
   G__memvar_setup((void*)((long)(&p->FindNum)-(long)(p)),105,0,0,-1,-1,-1,1,"FindNum=",0,"Finder Number");
   G__memvar_setup((void*)((long)(&p->ShowerHitNum)-(long)(p)),105,0,0,-1,-1,-1,1,"ShowerHitNum=",0,"Shower Hit Number");
   G__memvar_setup((void*)((long)(&p->ParCode)-(long)(p)),105,0,0,-1,-1,-1,1,"ParCode=",0,"Particle Code");
   G__memvar_setup((void*)((long)(&p->Fup)-(long)(p)),105,0,0,-1,-1,-1,1,"Fup=",0,"Finder set up  even 0 odd 1");
   G__memvar_setup((void*)((long)(&p->XX)-(long)(p)),102,0,0,-1,-1,-1,1,"XX=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->YY)-(long)(p)),102,0,0,-1,-1,-1,1,"YY=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->ZZ)-(long)(p)),102,0,0,-1,-1,-1,1,"ZZ=",0,"X, Y, Z positions in case of tracks , vertex in case of particle info");
   G__memvar_setup((void*)((long)(&p->pmag)-(long)(p)),102,0,0,-1,-1,-1,1,"pmag=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->pt)-(long)(p)),102,0,0,-1,-1,-1,1,"pt=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->pp)-(long)(p)),102,0,0,-1,-1,-1,1,"pp=",0,"in case of particle info, pmagnitude, ptheta and pphi");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__HitPosdictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarHitPosdict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncHitPos(void) {
   /* HitPos */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__HitPosdictLN_HitPos));
   G__memfunc_setup("HitPos",599,G__HitPosdict_233_0_1, 105, G__get_linked_tagnum(&G__HitPosdictLN_HitPos), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("HitPos",599,G__HitPosdict_233_0_2, 105, G__get_linked_tagnum(&G__HitPosdictLN_HitPos), -1, 0, 1, 1, 1, 0, "f - 'Float_t' 0 - random", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__HitPosdict_233_0_3, 85, G__get_linked_tagnum(&G__HitPosdictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&HitPos::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__HitPosdict_233_0_4, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&HitPos::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__HitPosdict_233_0_5, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&HitPos::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__HitPosdict_233_0_6, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&HitPos::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__HitPosdictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__HitPosdict_233_0_10, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__HitPosdict_233_0_11, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&HitPos::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__HitPosdict_233_0_12, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&HitPos::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__HitPosdict_233_0_13, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&HitPos::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__HitPosdict_233_0_14, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&HitPos::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("HitPos", 599, G__HitPosdict_233_0_15, (int) ('i'), G__get_linked_tagnum(&G__HitPosdictLN_HitPos), -1, 0, 1, 1, 1, 0, "u 'HitPos' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~HitPos", 725, G__HitPosdict_233_0_16, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__HitPosdict_233_0_17, (int) ('u'), G__get_linked_tagnum(&G__HitPosdictLN_HitPos), -1, 1, 1, 1, 1, 0, "u 'HitPos' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncHitPosdict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {
}

static void G__cpp_setup_global2() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalHitPosdict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
  G__cpp_setup_global2();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {
}

static void G__cpp_setup_func13() {
}

static void G__cpp_setup_func14() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcHitPosdict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
  G__cpp_setup_func13();
  G__cpp_setup_func14();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__HitPosdictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__HitPosdictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__HitPosdictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__HitPosdictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__HitPosdictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__HitPosdictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__HitPosdictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__HitPosdictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__HitPosdictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR = { "iterator<bidirectional_iterator_tag,TObject*,long,const TObject**,const TObject*&>" , 115 , -1 };
G__linked_taginfo G__HitPosdictLN_TVectorTlEfloatgR = { "TVectorT<float>" , 99 , -1 };
G__linked_taginfo G__HitPosdictLN_TVectorTlEdoublegR = { "TVectorT<double>" , 99 , -1 };
G__linked_taginfo G__HitPosdictLN_HitPos = { "HitPos" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableHitPosdict() {
  G__HitPosdictLN_TClass.tagnum = -1 ;
  G__HitPosdictLN_TBuffer.tagnum = -1 ;
  G__HitPosdictLN_TMemberInspector.tagnum = -1 ;
  G__HitPosdictLN_TObject.tagnum = -1 ;
  G__HitPosdictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__HitPosdictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__HitPosdictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__HitPosdictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__HitPosdictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR.tagnum = -1 ;
  G__HitPosdictLN_TVectorTlEfloatgR.tagnum = -1 ;
  G__HitPosdictLN_TVectorTlEdoublegR.tagnum = -1 ;
  G__HitPosdictLN_HitPos.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableHitPosdict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__HitPosdictLN_TClass);
   G__get_linked_tagnum_fwd(&G__HitPosdictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__HitPosdictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__HitPosdictLN_TObject);
   G__get_linked_tagnum_fwd(&G__HitPosdictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__HitPosdictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__HitPosdictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__HitPosdictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__HitPosdictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR);
   G__get_linked_tagnum_fwd(&G__HitPosdictLN_TVectorTlEfloatgR);
   G__get_linked_tagnum_fwd(&G__HitPosdictLN_TVectorTlEdoublegR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__HitPosdictLN_HitPos),sizeof(HitPos),-1,62720,"A track segment",G__setup_memvarHitPos,G__setup_memfuncHitPos);
}
extern "C" void G__cpp_setupHitPosdict(void) {
  G__check_setup_version(30051515,"G__cpp_setupHitPosdict()");
  G__set_cpp_environmentHitPosdict();
  G__cpp_setup_tagtableHitPosdict();

  G__cpp_setup_inheritanceHitPosdict();

  G__cpp_setup_typetableHitPosdict();

  G__cpp_setup_memvarHitPosdict();

  G__cpp_setup_memfuncHitPosdict();
  G__cpp_setup_globalHitPosdict();
  G__cpp_setup_funcHitPosdict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncHitPosdict();
  return;
}
class G__cpp_setup_initHitPosdict {
  public:
    G__cpp_setup_initHitPosdict() { G__add_setup_func("HitPosdict",(G__incsetup)(&G__cpp_setupHitPosdict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initHitPosdict() { G__remove_setup_func("HitPosdict"); }
};
G__cpp_setup_initHitPosdict G__cpp_setup_initializerHitPosdict;

