#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class HCircle++;
#pragma link C++ class HSprial++;
#pragma link C++ class HLine++;
#pragma link C++ class HTrack++;
#pragma link C++ class TPCHit++;
#pragma link C++ class TempCluster++;
#pragma link C++ class TPCCluster++;
#pragma link C++ class std::vector< HCircle >++;
#pragma link C++ class std::vector< HSprial >++;
#pragma link C++ class std::vector< HLine >++;
#pragma link C++ class std::vector< HCircle, allocator<HCircle> >::iterator++;
#pragma link C++ class std::vector< HSprial, allocator<HSprial> >::iterator++;
#pragma link C++ class std::vector< HLine, allocator<HLine> >::iterator++;
#pragma link C++ class std::vector< TPCHit >++;
#pragma link C++ class std::vector< TPCHit, allocator<TPCHit> >::iterator++;
#pragma link C++ class std::vector< TPCCluster >++;
#pragma link C++ class std::vector< TPCCluster, allocator<TPCCluster> >::iterator++;
#pragma link C++ class std::vector< std::vector<TPCCluster> >;

#pragma link C++ class std::vector< TPolyLine3D* >++;
#pragma link C++ class std::vector< TPolyLine3D*, allocator<TPolyLine3D* > >::iterator++;
#endif // __CINT__
