#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class  TPCPadHit+;
#pragma link C++ class  TPCPadHitCluster+;
//#pragma link C++ class  TPCTrack++;
#pragma link C++ function TPCClusterer;
#pragma link C++ function TPCFindBlock;
#pragma link C++ function TPCFindEdgeBlock;
#pragma link C++ function TPCFindBlockCircle;
#pragma link C++ function AdjacentCluster;
#pragma link C++ function ClusteringCalculator;
#pragma link C++ function DrawClusterArrXZ;
#pragma link C++ function DrawClusterArrYZ;
#pragma link C++ function CircleCrossing;
#pragma link C++ function yDistCalculator;
#pragma link C++ function TPCRootClusters;
#pragma link C++ function TPCEdgeClusters;
#pragma link C++ function TPCTrackBlocking;

#pragma link C++ class std::vector<TVector2>+;
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#pragma link C++ class std::vector< TPCPadHit >+;
#pragma link C++ class std::vector< TPCPadHitCluster >++;
#pragma link C++ class std::vector< TPCPadHitCluster, allocator<TPCPadHitCluster> >::iterator>+;
#pragma link C++ class std::vector< TPCPadHitCluster, allocator<TPCPadHitCluster> >::iterator>-;
#pragma link C++ class std::vector< TPCPadHitCluster, allocator<TPCPadHitCluster> >::iterator>++;
#pragma link C++ operators vector<vector<TPCPadHitCluster,allocator<TPCPadHitCluster> >::iterator>::iterator;
#pragma link C++ operators vector<vector<TPCPadHitCluster,allocator<TPCPadHitCluster> >::iterator>::const_iterator;
#pragma link C++ operators vector<vector<TPCPadHitCluster,allocator<TPCPadHitCluster> >::iterator>::reverse_iterator;

#endif // __CINT__
