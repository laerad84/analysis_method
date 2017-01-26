void Test1(){

  TClonesArray* arr = new TClonesArray("TGraph");
  for( int i = 0; i< 10; i++){
    new((*arr)[i]) TGraph();
    TGraph* gr = (TGraph*)arr->At(i);
    gr->SetUniqueID(i);
  }
  arr->RemoveAt(5);
  std::cout<< arr->GetEntries() << std::endl;
  for( int i = 0; i< arr->GetEntries(); i++){
    TGraph* gr = (TGraph*)arr->At(i);
    if( gr == NULL ){ std::cout<< i << "NULL" << std::endl; continue; }
    std::cout<< ((TGraph*)arr->At(i))->GetUniqueID() << std::endl;

  }
}
