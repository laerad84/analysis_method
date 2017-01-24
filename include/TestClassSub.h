#ifndef TestClassSub__H__
#define TestClassSub__H__
#include "TestClass.h"
#include "TObject.h"

class TestClassSub {
 public:
  TestClassSub();
  ~TestClassSub();

  TestClass* cl;
  void Run();
  ClassDef(TestClassSub, 0 )
};

#endif
