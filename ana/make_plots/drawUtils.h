void applyStyle(TGraph* myhist){

  myhist->GetXaxis()->CenterTitle(true);
  myhist->GetYaxis()->CenterTitle(true);

  myhist->SetMarkerStyle(8);
  myhist->SetLineWidth(4);//(2);

  //  myhist->SetLabelFont(42);
  //  myhist->SetTitleFont(42);

  //myhist->GetZaxis()->SetTitleFont(42);
  myhist->GetYaxis()->SetTitleFont(42);
  myhist->GetXaxis()->SetTitleFont(42);

  myhist->GetXaxis()->SetLabelSize(0.05);
  myhist->GetXaxis()->SetTitleSize(0.05);
  myhist->GetYaxis()->SetLabelSize(0.05);
  myhist->GetYaxis()->SetTitleSize(0.05);
  //myhist->GetZaxis()->SetLabelSize(0.05);
  // myhist->GetZaxis()->SetTitleSize(0.05);
  
  myhist->GetXaxis()->SetTickLength(0.01);
  myhist->GetYaxis()->SetTickLength(0.02);
  
  myhist->GetXaxis()->SetNdivisions(510);
  myhist->GetYaxis()->SetNdivisions(510);
}


void applyStyle(MnvH1D* myhist){

  myhist->GetXaxis()->CenterTitle(true);
  myhist->GetYaxis()->CenterTitle(true);

  myhist->SetMarkerStyle(8);
  myhist->SetLineWidth(4);//(2);

  myhist->SetLabelFont(42);
  myhist->SetTitleFont(42);

  myhist->GetYaxis()->SetTitleFont(42);
  myhist->GetXaxis()->SetTitleFont(42);

  myhist->SetLabelSize(0.05,"x");
  myhist->SetTitleSize(0.05,"x");
  myhist->SetLabelSize(0.05,"y");
  myhist->SetTitleSize(0.05,"y");
  
  myhist->SetTickLength(0.01, "Y");
  myhist->SetTickLength(0.02, "X");
  
  myhist->SetNdivisions(510, "XYZ");
}

void applyStyle(MnvH2D* myhist){

  myhist->GetXaxis()->CenterTitle(true);
  myhist->GetYaxis()->CenterTitle(true);
  myhist->GetZaxis()->CenterTitle(true);
  myhist->SetMarkerStyle(8);
  myhist->SetLineWidth(4);//(2);

  myhist->SetLabelFont(42);
  myhist->SetTitleFont(42);

  myhist->GetYaxis()->SetTitleFont(42);
  myhist->GetXaxis()->SetTitleFont(42);

  
  myhist->SetLabelSize(0.05,"x");
  myhist->SetTitleSize(0.05,"x");
  myhist->SetLabelSize(0.05,"y");
  myhist->SetTitleSize(0.05,"y");  
  myhist->SetLabelSize(0.05,"z");
  myhist->SetTitleSize(0.04,"z");
  myhist->GetZaxis()->SetTitleOffset(1.05);

  myhist->SetTickLength(0.01, "Y");
  myhist->SetTickLength(0.02, "X");
  
  myhist->SetNdivisions(510, "XYZ");
}

void applyStyle(TH1D* myhist){

  myhist->GetXaxis()->CenterTitle(true);
  myhist->GetYaxis()->CenterTitle(true);

  myhist->SetMarkerStyle(8);
  myhist->SetLineWidth(4);//(2);

  myhist->SetLabelFont(42);
  myhist->SetTitleFont(42);

  myhist->GetYaxis()->SetTitleFont(42);
  myhist->GetXaxis()->SetTitleFont(42);

  myhist->SetLabelSize(0.05,"x");
  myhist->SetTitleSize(0.05,"x");
  myhist->SetLabelSize(0.05,"y");
  myhist->SetTitleSize(0.05,"y");
  
  myhist->SetTickLength(0.01, "Y");
  myhist->SetTickLength(0.02, "X");
  
  myhist->SetNdivisions(510, "XYZ");
}

void drawVariableQELike(string varname, TFile *myfile,bool dolog, string xtitle, string ytitle, string title,double yscale){
  TCanvas *c1 = new TCanvas(varname.c_str(),varname.c_str());
  MnvH1D *qe = (MnvH1D*)myfile->Get(Form("%s_qelike_qe",varname.c_str()));
  MnvH1D *res = (MnvH1D*)myfile->Get(Form("%s_qelike_res",varname.c_str()));
  MnvH1D *dis = (MnvH1D*)myfile->Get(Form("%s_qelike_dis",varname.c_str()));
  MnvH1D *npnh = (MnvH1D*)myfile->Get(Form("%s_qelike_oth",varname.c_str()));
  MnvH1D *bkg = (MnvH1D*)myfile->Get(Form("%s_qelikenot",varname.c_str()));
  MnvH1D *data = (MnvH1D*)myfile->Get(Form("%s_data",varname.c_str()));

  qe->SetFillColor(kMagenta-4);
  res->SetFillColor(kPink+9);
  dis->SetFillColor(kMagenta-9);
  npnh->SetFillColor(kBlack);
  bkg->SetFillColor(kBlue-7);

  THStack *hs = new THStack();
  hs->Add(new TH1D(bkg->GetCVHistoWithError()));
  hs->Add(new TH1D(dis->GetCVHistoWithError()));
  hs->Add(new TH1D(res->GetCVHistoWithError()));
  hs->Add(new TH1D(npnh->GetCVHistoWithError()));
  hs->Add(new TH1D(qe->GetCVHistoWithError()));

  TLegend *leg = new TLegend(0.55,0.7,0.85,0.91);
  leg->AddEntry(data,"MINERvA Data","p");
  leg->AddEntry(qe,"QELike-QE","F");
  leg->AddEntry(npnh,"QELike-2p2h","F");
  leg->AddEntry(res,"QELike-Res","F");
  leg->AddEntry(dis,"QELike-DIS","F");
  leg->AddEntry(bkg,"Background","F");

  data->GetXaxis()->SetTitle(xtitle.c_str());
  data->GetYaxis()->SetTitle(ytitle.c_str());
  data->SetTitle(title.c_str());
  applyStyle(data);

  data->Draw("PE1");
  hs->Draw("HISTSAME");
  data->Draw("PE1SAME");
  dolog?data->GetYaxis()->SetRangeUser(0.1,yscale*data->GetBinContent(data->GetMaximumBin())):data->GetYaxis()->SetRangeUser(0,yscale*data->GetBinContent(data->GetMaximumBin()));
  leg->Draw("SAME");
  if(dolog) c1->SetLogy(true);
  gPad->RedrawAxis();
  dolog?c1->Print(Form("%s_log.png",varname.c_str())):c1->Print(Form("%s.png",varname.c_str()));
  dolog?c1->Print(Form("%s_log.C",varname.c_str())):c1->Print(Form("%s.C",varname.c_str()));
  dolog?c1->Print(Form("%s_log.eps",varname.c_str())):c1->Print(Form("%s.eps",varname.c_str()));

}


void drawVariableQELikeProjX(string varname, TFile *myfile,bool dolog, string xtitle, string ytitle, string title,double yscale){
  TCanvas *c1 = new TCanvas(varname.c_str(),varname.c_str());
  MnvH1D *qe = (MnvH1D*)((MnvH2D*)myfile->Get(Form("%s_qelike_qe",varname.c_str())))->ProjectionX();
  MnvH1D *res = (MnvH1D*)((MnvH2D*)myfile->Get(Form("%s_qelike_res",varname.c_str())))->ProjectionX();
  MnvH1D *dis = (MnvH1D*)((MnvH2D*)myfile->Get(Form("%s_qelike_dis",varname.c_str())))->ProjectionX();
  MnvH1D *npnh = (MnvH1D*)((MnvH2D*)myfile->Get(Form("%s_qelike_oth",varname.c_str())))->ProjectionX();
  MnvH1D *bkg = (MnvH1D*)((MnvH2D*)myfile->Get(Form("%s_qelikenot",varname.c_str())))->ProjectionX();
  MnvH1D *data = (MnvH1D*)((MnvH2D*)myfile->Get(Form("%s_data",varname.c_str())))->ProjectionX();
  
  qe->SetFillColor(kMagenta-4);
  res->SetFillColor(kGreen-3);
  dis->SetFillColor(kOrange-2);
  npnh->SetFillColor(kViolet-5);
  bkg->SetFillColor(kBlue-7);


  qe->SetLineColor(kMagenta-4);
  res->SetLineColor(kGreen-3);
  dis->SetLineColor(kOrange-2);
  npnh->SetLineColor(kViolet-5);
  bkg->SetLineColor(kBlue-7);

  
  THStack *hs = new THStack();
  hs->Add(new TH1D(bkg->GetCVHistoWithError()));
  hs->Add(new TH1D(npnh->GetCVHistoWithError()));
  hs->Add(new TH1D(dis->GetCVHistoWithError()));
  hs->Add(new TH1D(res->GetCVHistoWithError()));
  hs->Add(new TH1D(qe->GetCVHistoWithError()));

  TLegend *leg = new TLegend(0.55,0.7,0.85,0.91);
  leg->SetFillColor(kWhite);
  leg->AddEntry(data,"MINERvA Data","p");
  leg->AddEntry(qe,"QELike-QE","F");
  leg->AddEntry(npnh,"QELike-2p2h","F");
  leg->AddEntry(res,"QELike-Res","F");
  leg->AddEntry(dis,"QELike-DIS","F");
  leg->AddEntry(bkg,"Background","F");

  data->GetXaxis()->SetTitle(xtitle.c_str());
  data->GetYaxis()->SetTitle(ytitle.c_str());
  data->SetTitle(title.c_str());
  applyStyle(data);


  data->GetXaxis()->SetNoExponent(true);
  data->GetXaxis()->SetNdivisions(505);


  data->Draw("PE1");
  hs->Draw("HISTSAME");
  data->Draw("PE1SAME");
  dolog?data->GetYaxis()->SetRangeUser(0.1,yscale*data->GetBinContent(data->GetMaximumBin())):data->GetYaxis()->SetRangeUser(0,yscale*data->GetBinContent(data->GetMaximumBin()));
  leg->Draw("SAME");
  if(dolog) c1->SetLogy(true);
  gPad->RedrawAxis();
  dolog?c1->Print(Form("%s_log.png",varname.c_str())):c1->Print(Form("%s.png",varname.c_str()));
  dolog?c1->Print(Form("%s_log.C",varname.c_str())):c1->Print(Form("%s.C",varname.c_str()));
  dolog?c1->Print(Form("%s_log.eps",varname.c_str())):c1->Print(Form("%s.eps",varname.c_str()));

}

void drawVariableFSPion(string varname, TFile *myfile,bool dolog, string xtitle, string ytitle, string title,double yscale,string specname = "NULL",bool binwidthnorm = false){
  TCanvas *c1 = new TCanvas(varname.c_str(),varname.c_str());
  MnvH1D *qe = (MnvH1D*)myfile->Get(Form("%s_qelike",varname.c_str()));
  MnvH1D *sch = (MnvH1D*)myfile->Get(Form("%s_qelikenot_singlechargedpion",varname.c_str()));
  MnvH1D *snu = (MnvH1D*)myfile->Get(Form("%s_qelikenot_singleneutralpion",varname.c_str()));
  MnvH1D *multi = (MnvH1D*)myfile->Get(Form("%s_qelikenot_multipion",varname.c_str()));
  MnvH1D *other = (MnvH1D*)myfile->Get(Form("%s_qelikenot_nopions",varname.c_str()));
  MnvH1D *data = (MnvH1D*)myfile->Get(Form("%s_data",varname.c_str()));
  qe->SetFillColor(kMagenta-4);
  sch->SetFillColor(kBlue-7);
  snu->SetFillColor(kGreen-3);
  multi->SetFillColor(kRed-4);
  other->SetFillColor(kBlack);

  qe->SetLineColor(kMagenta-4);
  sch->SetLineColor(kBlue-7);
  snu->SetLineColor(kGreen-3);
  multi->SetLineColor(kRed-4);
  other->SetLineColor(kBlack);

  
  if(binwidthnorm){
    qe->Scale(1.0,"width");
    sch->Scale(1.0,"width");
    snu->Scale(1.0,"width");
    multi->Scale(1.0,"width");
    other->Scale(1.0,"width");
    data->Scale(1.0,"width");
  }

  THStack *hs = new THStack();
  hs->Add(new TH1D(other->GetCVHistoWithError()));
  hs->Add(new TH1D(multi->GetCVHistoWithError()));
  hs->Add(new TH1D(snu->GetCVHistoWithError()));
  hs->Add(new TH1D(sch->GetCVHistoWithError()));
  hs->Add(new TH1D(qe->GetCVHistoWithError()));

  TLegend *leg = new TLegend(0.55,0.7,0.85,0.91);
  leg->SetFillColor(kWhite);
  leg->AddEntry(data,"MINERvA Data","p");
  leg->AddEntry(qe,"QELike","F");
  leg->AddEntry(sch,"Single #pi^{+/-} in FS","F");
  leg->AddEntry(snu,"Single #pi^{0} in FS","F");
  leg->AddEntry(multi,"N#pi in FS","F");
  leg->AddEntry(other,"Other","F");


  data->GetXaxis()->SetTitle(xtitle.c_str());
  data->GetYaxis()->SetTitle(ytitle.c_str());
  data->GetYaxis()->SetTitleOffset(1.4);
  data->SetTitle(title.c_str());
  applyStyle(data);

  data->Draw("PE1");
  hs->Draw("HISTSAME");
  data->Draw("PE1SAME");
  dolog?data->GetYaxis()->SetRangeUser(0.1,yscale*data->GetBinContent(data->GetMaximumBin())):data->GetYaxis()->SetRangeUser(0,yscale*data->GetBinContent(data->GetMaximumBin()));
  leg->Draw("SAME");
  if(dolog) c1->SetLogy(true);
  gPad->RedrawAxis();
  if(specname!="NULL") varname = varname+specname;
  dolog?c1->Print(Form("%s_log.png",varname.c_str())):c1->Print(Form("%s.png",varname.c_str()));
  dolog?c1->Print(Form("%s_log.C",varname.c_str())):c1->Print(Form("%s.C",varname.c_str()));
  dolog?c1->Print(Form("%s_log.eps",varname.c_str())):c1->Print(Form("%s.eps",varname.c_str()));

}



