void Test(){

  gSystem->Load("lib/libanalysis_method.dylib");


  HSprial sprial;
  sprial.R = 2;
  sprial.DY = 2;
  sprial.DTheta = 2*TMath::Pi();
  sprial.X=0;
  sprial.Z=0;
  sprial.RL=1;
  sprial.InitPos=TVector3(2,0,0);
  TPolyLine3D* line = sprial.GenerateHelix();
  line->Draw();
  Double_t dR,dY;
  sprial.CalculateDist(TVector3(-2.1,0.5,0),dR,dY);
  std::cout<< dR << "\t" << dY << std::endl;
}
