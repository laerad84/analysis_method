void testVector(){

  std::vector<int> vec;
  for( int i = 0; i< 10; i++){
    vec.push_back(i);
  }

  vec.erase(vec.begin()+5);
  for( int i = 0; i< vec.size(); i++){
    std::cout<< vec.at(i) << "\t";
  }
  std::cout<< std::endl;

  std::vector<int> vec1;
  vec1.push_back(3);
  vec.erase(vec.begin()+3);

  for( int ivec = 0; ivec < vec1.size(); ivec++){
    for( int iivec = 0; iivec< vec.size(); ){
      if(TMath::Abs(vec1.at( ivec ) - vec.at(iivec) )  < 2 ){
	vec1.push_back( vec.at(iivec));
	vec.erase( vec.begin() + iivec );
      }else{
	iivec++;
      }
    }
  }

  for( int ivec = 0; ivec < vec1.size(); ivec++){
    std::cout<< ivec << "\t" << vec1.at(ivec) << std::endl;
  }

}
