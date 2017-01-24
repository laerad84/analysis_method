#include <iostream>

class TestClass {
 public:
  TestClass(){ std::cout<< "test" << std::endl;}
  ~TestClass(){ std::cout<< "distruct" << std::endl; }

  void Run() { std::cout<< "Run" << std::endl;}
};
