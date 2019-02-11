// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME srcdICParticle_dict

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
#include "inc/CParticle.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_CParticle(void *p = 0);
   static void *newArray_CParticle(Long_t size, void *p);
   static void delete_CParticle(void *p);
   static void deleteArray_CParticle(void *p);
   static void destruct_CParticle(void *p);
   static void streamer_CParticle(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CParticle*)
   {
      ::CParticle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::CParticle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("CParticle", ::CParticle::Class_Version(), "inc/CParticle.h", 17,
                  typeid(::CParticle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::CParticle::Dictionary, isa_proxy, 16,
                  sizeof(::CParticle) );
      instance.SetNew(&new_CParticle);
      instance.SetNewArray(&newArray_CParticle);
      instance.SetDelete(&delete_CParticle);
      instance.SetDeleteArray(&deleteArray_CParticle);
      instance.SetDestructor(&destruct_CParticle);
      instance.SetStreamerFunc(&streamer_CParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CParticle*)
   {
      return GenerateInitInstanceLocal((::CParticle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::CParticle*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr CParticle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *CParticle::Class_Name()
{
   return "CParticle";
}

//______________________________________________________________________________
const char *CParticle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CParticle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int CParticle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CParticle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *CParticle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CParticle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *CParticle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CParticle*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void CParticle::Streamer(TBuffer &R__b)
{
   // Stream an object of class CParticle.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> Px;
      R__b >> Py;
      R__b >> Pz;
      R__b >> E;
      R__b >> Mass;
      R__b >> Charge;
      R__b >> m_pid;
      R__b >> m_status;
      {
         vector<int> &R__stl =  parameters;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            int R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      R__b.CheckByteCount(R__s, R__c, CParticle::IsA());
   } else {
      R__c = R__b.WriteVersion(CParticle::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << Px;
      R__b << Py;
      R__b << Pz;
      R__b << E;
      R__b << Mass;
      R__b << Charge;
      R__b << m_pid;
      R__b << m_status;
      {
         vector<int> &R__stl =  parameters;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<int>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_CParticle(void *p) {
      return  p ? new(p) ::CParticle : new ::CParticle;
   }
   static void *newArray_CParticle(Long_t nElements, void *p) {
      return p ? new(p) ::CParticle[nElements] : new ::CParticle[nElements];
   }
   // Wrapper around operator delete
   static void delete_CParticle(void *p) {
      delete ((::CParticle*)p);
   }
   static void deleteArray_CParticle(void *p) {
      delete [] ((::CParticle*)p);
   }
   static void destruct_CParticle(void *p) {
      typedef ::CParticle current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_CParticle(TBuffer &buf, void *obj) {
      ((::CParticle*)obj)->::CParticle::Streamer(buf);
   }
} // end of namespace ROOT for class ::CParticle

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 210,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace {
  void TriggerDictionaryInitialization_CParticle_dict_Impl() {
    static const char* headers[] = {
"inc/CParticle.h",
0
    };
    static const char* includePaths[] = {
"/users/admin/share/sl7/root_v6.12.06_gcc485/include",
"/users/whopkins/FastSimNN/validate/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "CParticle_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Integrate this class into ROOT)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Integrate this class into ROOT)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Integrate this class into ROOT)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Integrate this class into ROOT)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$inc/CParticle.h")))  CParticle;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "CParticle_dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "inc/CParticle.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"CParticle", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("CParticle_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_CParticle_dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_CParticle_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_CParticle_dict() {
  TriggerDictionaryInitialization_CParticle_dict_Impl();
}
