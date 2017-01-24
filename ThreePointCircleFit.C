void ThreePointCircleFit(){

  TVector2 vec[3];

  for( int i = 0; i< 3; i++){
    double x,y;
    gRandom->Rannor(x,y);
    vec[i]=TVector2(x,y);
    std::cout<< x << "\t" << y << "\t" << vec[i].X() << "\t" << vec[i].Y() << std::endl;
  }

  Double_t meanX[3];
  Double_t meanY[3];
  Double_t deltaX[3];
  Double_t deltaY[3];
  std::cout<< "Mean and delta " << std::endl;
  for( int i = 0; i< 3; i++){
    std::cout<< (i+1)%3 << std::endl;
    meanX[i] = (vec[i].X() + vec[(i+1)%3].X())*0.5;
    meanY[i] = (vec[i].Y() + vec[(i+1)%3].Y())*0.5;
    deltaX[i] = (vec[(i+1)%3].X() - vec[i].X());
    deltaY[i] = (vec[(i+1)%3].Y() - vec[i].Y());
    std::cout<< meanX[i] << "\t" << meanY[i] << "\t" << deltaX[i] << "\t" << deltaY[i] << std::endl;
  }

  Double_t centerX(0),centerY(0);
  centerX =
    (meanX[1]*deltaX[1]*deltaY[0] - meanX[0]*deltaX[0]*deltaY[1] - 0.5*deltaY[0]*deltaY[1]*deltaY[2])/
    (deltaY[0]*deltaX[1] - deltaY[1]*deltaX[0]);
  centerY =
    (meanY[1]*deltaY[1]*deltaX[0] - meanY[0]*deltaY[0]*deltaX[1] - 0.5*deltaX[0]*deltaX[1]*deltaX[2])/
    (deltaX[0]*deltaY[1] - deltaX[1]*deltaY[0]);
  std::cout<< centerX << "\t" << centerY << std::endl;


  Double_t r[3];
  for( int i = 0; i< 3; i++){
    r[i] = TMath::Sqrt( TMath::Power(vec[i].X() - centerX ,2 ) + TMath::Power(vec[i].Y()-centerY,2));
    std::cout << r[i] << std::endl;
  }

  TGraph* gr = new TGraph();
  gr->SetMarkerStyle(22);
  TMarker* marker = new TMarker(centerX, centerY, 21);
  for( int i = 0; i< 3; i++){
    gr->SetPoint( i, vec[i].X(), vec[i].Y());
  }



  TCanvas* can = new TCanvas("can","can",800,800);
  can->DrawFrame( centerX+r[0]*-1.1, centerY+r[0]*-1.1,centerX+r[0]*1.1,centerY+r[0]*1.1);

  TArc* arc = new TArc(centerX,centerY,r[0]);
  arc->SetFillColor(0);
  arc->SetFillStyle(0);
  arc->Draw();
  gr->Draw("p");
  marker->Draw("p");

}
