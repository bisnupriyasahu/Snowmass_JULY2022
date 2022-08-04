{
  auto *c1 = new TCanvas("c1","",800,800);
  THStack *stackHisto = new THStack("stackHisto","Title");
  TFile* datafile1 = TFile::Open("TTBAR.root");
  TFile* datafile2 = TFile::Open("DY.root");
  //  TFile* datafile3 = TFile::Open("case2_Had_output_300_1_L.root");
  //TFile* datafile4 = TFile::Open("case2_Had_output_300_1_R.root");

  //TFile* datafile5 = TFile::Open("case2_Had_output_500_350_L.root");
  //TFile* datafile6 = TFile::Open("case2_Had_output_500_350_R.root");

  TFile* datafile7 = TFile::Open("case2_output_1000_1_L.root");
  TFile* datafile8 = TFile::Open("case2_output_1000_1_R.root");

  TFile* datafile9 = TFile::Open("case2_output_800_300_L.root");
  TFile* datafile10 = TFile::Open("case2_output_800_300_R.root");

  TH1F *h1 = dynamic_cast<TH1F*>(datafile1->Get("MT"));
  TH1F *h2 = dynamic_cast<TH1F*>(datafile2->Get("MT"));
  //TH1F *h3 = dynamic_cast<TH1F*>(datafile3->Get("MT"));
  //TH1F *h4 = dynamic_cast<TH1F*>(datafile4->Get("MT"));

  
  // TH1F *h5 = dynamic_cast<TH1F*>(datafile5->Get("MT"));
  //TH1F *h6 = dynamic_cast<TH1F*>(datafile6->Get("MT"));
  TH1F *h7 = dynamic_cast<TH1F*>(datafile7->Get("MT"));
  TH1F *h8 = dynamic_cast<TH1F*>(datafile8->Get("MT"));

  TH1F *h9 = dynamic_cast<TH1F*>(datafile9->Get("MT"));
  TH1F *h10 = dynamic_cast<TH1F*>(datafile10->Get("MT"));

  TH1F *h_n1 = dynamic_cast<TH1F*>(datafile1->Get("nEvents"));
  TH1F *h_n2 = dynamic_cast<TH1F*>(datafile2->Get("nEvents"));
  //TH1F *h_n3 = dynamic_cast<TH1F*>(datafile3->Get("nEvents_1"));
  //TH1F *h_n4 = dynamic_cast<TH1F*>(datafile4->Get("nEvents_1"));

  
  //TH1F *h_n5 = dynamic_cast<TH1F*>(datafile5->Get("nEvents_1"));
  //TH1F *h_n6 = dynamic_cast<TH1F*>(datafile6->Get("nEvents_1"));
  TH1F *h_n7 = dynamic_cast<TH1F*>(datafile7->Get("nEvents_1"));
  TH1F *h_n8 = dynamic_cast<TH1F*>(datafile8->Get("nEvents_1"));

  TH1F *h_n9 = dynamic_cast<TH1F*>(datafile9->Get("nEvents_1"));
  TH1F *h_n10 = dynamic_cast<TH1F*>(datafile10->Get("nEvents_1"));


  //  TH1F *h6 = dynamic_cast<TH1F*>(datafile3->Get("nEvents"));

  float luminosity = 3000000; // 1/pb
  float TTXsec = 864.4; //pb
  float DYXsec = 246.5; //pb
  float sXsec_300 =105.081; //pb
  float sXsec_500 = 6.68102; //pb
  float sXsec_800 = 0.387393; //pb
  float sXsec_1000 = 0.0875693; //pb



  TTEvts = h_n1->GetBinContent(2);
  DYEvts = h_n2->GetBinContent(2);
    
  sEvts_1000_1_L = h_n7->GetBinContent(2);
  sEvts_1000_1_R = h_n8->GetBinContent(2);
  
  sEvts_800_300_L = h_n9->GetBinContent(2);
  sEvts_800_300_R = h_n10->GetBinContent(2);

  h1->SetFillColor(kGreen);
  h2->SetFillColor(kBlue);

  
  h7->SetLineColor(kRed);
  h8->SetLineColor(kYellow);

  h9->SetLineColor(kBlack);
  h10->SetLineColor(kCyan);

  //  h7->SetLineStyle(10);                                                                                                                                                            
  // h8->SetLineStyle(3);                                                                                                                                                              
  //h9->SetLineStyle(10);                                                                                                                                                              
  h10->SetLineStyle(9);


  
  
  h7->SetLineWidth(2);
  h8->SetLineWidth(2);
  h9->SetLineWidth(2);
  h10->SetLineWidth(2);

 
  

  std::cout<<"weight of TTBAR : "<< luminosity * TTXsec / TTEvts<<std::endl;
  std::cout<<"weight of DY : "<< luminosity * DYXsec / DYEvts<< "  "<<DYXsec<< " "<< DYEvts <<std::endl;

  std::cout<<"h1 Integral before : "<< h1->Integral()<<std::endl;
  std::cout<<"h2 Integral before : "<< h2->Integral()<<std::endl;
   std::cout<<"h7 Integral before : "<< h7->Integral()<<std::endl;
  std::cout<<"h8 Integral before : "<< h8->Integral()<<std::endl;
  std::cout<<"h9 Integral before : "<< h9->Integral()<<std::endl;
  std::cout<<"h10 Integral before : "<< h10->Integral()<<std::endl;
  std::cout<<"-------------------------------------------------------------------"<<std::endl;
  
    h1->Scale(luminosity * TTXsec / TTEvts);
    h2->Scale(luminosity * DYXsec / DYEvts);
  
    h7->Scale(luminosity * sXsec_1000 / sEvts_1000_1_L );
    h8->Scale(luminosity * sXsec_1000 / sEvts_1000_1_R );

    h9->Scale(luminosity * sXsec_800 / sEvts_800_300_L );
    h10->Scale(luminosity * sXsec_800 / sEvts_800_300_R );

  std::cout<<"h1 Integral after : "<< h1->Integral()<<std::endl;
  std::cout<<"h2 Integral after : "<< h2->Integral()<<std::endl;
  std::cout<<"h7 Integral after : "<< h7->Integral()<<std::endl;
  std::cout<<"h8 Integral after : "<< h8->Integral()<<std::endl;
  std::cout<<"h9 Integral after : "<< h9->Integral()<<std::endl;
  std::cout<<"h10 Integral after : "<< h10->Integral()<<std::endl;

  float S1 = h7->Integral();
  float S2 = h8->Integral();
  float S3 = h9->Integral();
  float S4 = h10->Integral();
  float B = h1->Integral()+h2->Integral();

  float sig_1000_1_L = S1/sqrt(S1+B);
  float sig_1000_1_R = S2/sqrt(S2+B);
  float sig_800_300_L = S3/sqrt(S3+B);
  float sig_800_300_R = S4/sqrt(S4+B);
  std::cout<<"-------------------------------------------------------------------"<<std::endl;
  std::cout<<"S/sqrt(S+B) for left handed signal, mass point [STOP,LSP]: [1000,1]  : "<< sig_1000_1_L <<std::endl;
  std::cout<<"S/sqrt(S+B) for right handed signal, mass point [STOP,LSP]: [1000,1]  : "<< sig_1000_1_R <<std::endl;
  std::cout<<"S/sqrt(S+B) for left handed signal, mass point [STOP,LSP]: [800,300]  : "<< sig_800_300_L <<std::endl;
  std::cout<<"S/sqrt(S+B) for right handed signal, mass point [STOP,LSP]: [800,300]  : "<< sig_800_300_R <<std::endl;

    /* std::cout<<"h1 afterscaling : "<< h1->GetEntries()<<std::endl;
  std::cout<<"h2 afterscaling : "<< h2->GetEntries()<<std::endl;
  std::cout<<"h3 afterscaling : "<< h3->GetEntries()<<std::endl;
  std::cout<<"h4 afterscaling : "<< h4->GetEntries()<<std::endl;
    */
    //
    ///////////////////////////////////////////
    // h1->Draw();
    /*h1->Rebin(8);
    h2->Rebin(8);
    h7->Rebin(8);
    h8->Rebin(8);
    h9->Rebin(8);
   h10->Rebin(8);
    *///    stackHisto->Add(h3);  
  //stackHisto->Add(h4);
  stackHisto->Add(h2);
  stackHisto->Add(h1);

  std::cout<<"coming before : "<<std::endl;

  //h6->Draw("HIST");  
  //
  stackHisto->SetTitle("MT - Hadronically decayed #tau");  
  stackHisto->Draw("HIST");
  h7->Draw("HIST SAME");  
  h8->Draw("HIST SAME"); 
     
  h9->Draw("HIST SAME");
  h10->Draw("HIST SAME");

  
  // h5->Draw("HIST SAME");  
  //h6->Draw("HIST SAME");
  //h8->Draw("HIST SAME");  

  //  h9->Draw("HIST SAME");
  //h10->Draw("HIST SAME");  
  
  //    stackHisto->Draw("HIST");
  
  //stackHisto->GetXaxis()->SetTitle("Date [month/year]");
   std::cout<<"coming before 1A : "<<std::endl;    
  //  stackHisto->SetStats(0);
  //stackHisto->SetTitle("");
  stackHisto->GetXaxis()->SetRangeUser(0, 800);
  stackHisto->GetXaxis()->SetTitle("MT");
 std::cout<<"coming before 2A : "<<std::endl;    
  stackHisto->GetYaxis()->SetTitle("Events(GeV)");
    //stackHisto->GetXaxis()->SetLabelFont(42);
    //    stackHisto->GetXaxis()->SetLabelSize(0.06);
  

  stackHisto->GetXaxis()->SetTitleFont(42);
    stackHisto->GetXaxis()->SetTitleSize(0.045);

    //stackHisto->GetYaxis()->SetLabelFont(42);
    // stackHisto->GetYaxis()->SetLabelSize(0.06);
    stackHisto->GetYaxis()->SetTitleFont(42);
    stackHisto->GetYaxis()->SetTitleSize(0.045);
    stackHisto->GetYaxis()->SetTitleOffset(0.9);
    h2->SetLineColor(kBlack);
    h1->SetLineColor(kBlack);
       std::cout<<"coming before 2A : "<<std::endl;

  //

  // h1->SetLineWidth(2);
  // h2->SetLineWidth(2);

   //   h2->GetZaxis()->SetLabelFont(42);
   //h2->GetZaxis()->SetLabelSize(0.035);
   //h2->GetZaxis()->SetTitleSize(0.035);
   //h2->GetZaxis()->SetTitleFont(42);
 
       auto leg = new TLegend(0.6,0.6,0.8,0.9);
       //auto leg = new TLegend(0.1,0.7,0.6,0.9);
       //leg->AddEntry(eg_def,"default Ecaliso","l");
  leg->AddEntry(h1,"TTBar","F");
  leg->AddEntry(h2,"DY","F");
  /* leg->AddEntry(h3,"Signal[300,1] P_{#tau} = -1","l");
  leg->AddEntry(h4,"Signal[300,1] P_{#tau} = +1","l");
  leg->AddEntry(h5,"Signal[500,350] P_{#tau} = -1","l");
  leg->AddEntry(h6,"Signal[500,350] P_{#tau} = +1","l");
  */leg->AddEntry(h7,"Signal[1000,1] P_{#tau} = -1","l");
  leg->AddEntry(h8,"Signal[1000,1] P_{#tau} = +1","l");
   leg->AddEntry(h9,"Signal[800,300] P_{#tau} = -1","l");
  leg->AddEntry(h10,"Signal[800,300] P_{#tau} = +1","l");
  
  
  leg->SetFillColor(kWhite);
  leg->SetShadowColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.0300);
  leg->Draw();

  //  stackHisto->Draw();    
  c1->Modified();
    c1->SetLogy(1);
  c1->SaveAs("HAD_MT.png");
 
}
