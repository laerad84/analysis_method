#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class HCircle++;
#pragma link C++ class HHelix++;
#pragma link C++ class HLine++;
#pragma link C++ class HTrack++;
#pragma link C++ class TPCHit++;
#pragma link C++ class TPCCluster++;
#pragma link C++ class TPCTrack++;

#pragma link C++ class HPiM++;
#pragma link C++ class HKaonP++;
#pragma link C++ class HProton++;
#pragma link C++ class HLambda++;
#pragma link C++ class HCascade++;
#pragma link C++ class HDibaryonLL++;
#pragma link C++ class HDibaryonPC++;

#pragma link C++ class std::vector< HCircle >++;
#pragma link C++ class std::vector< HHelix >++;
#pragma link C++ class std::vector< HLine >++;
#pragma link C++ class std::vector< HCircle, allocator<HCircle> >::iterator++;
#pragma link C++ class std::vector< HHelix, allocator<HHelix> >::iterator++;
#pragma link C++ class std::vector< HLine, allocator<HLine> >::iterator++;
#pragma link C++ class std::vector< TPCHit >++;
#pragma link C++ class std::vector< TPCHit, allocator<TPCHit> >::iterator++;

#pragma link C++ class std::vector< TPCCluster >++;
#pragma link C++ class std::vector< TPCCluster, allocator<TPCCluster> >::iterator++;
#pragma link C++ class std::vector< std::vector<TPCCluster> >;

#pragma link C++ class std::vector< TPCTrack >++;

#pragma link C++ class std::vector< TPolyLine3D* >++;
#pragma link C++ class std::vector< TPolyLine3D*, allocator<TPolyLine3D* > >::iterator++;
#endif // __CINT__
