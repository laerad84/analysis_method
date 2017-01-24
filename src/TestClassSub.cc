#include "TestClassSub.h"
ClassImp(TestClassSub)
TestClassSub::TestClassSub(){
  cl = new TestClass();
}
TestClassSub::~TestClassSub(){
  delete cl;
}
void TestClassSub::Run(){
  cl->Run();
}
