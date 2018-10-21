void eventselectionplot(string htype="hx",float ymax=800.0, float xlmin=0.35, float ylmin=0.58){

string poln="LR";

	//TFile *fout = new TFile("eventselection.root","RECREATE");
TCanvas *c1 = new TCanvas("c1","multipads",800,600);
//TFile *fLR = new TFile( ("out_"+poln+".root").c_str());
//TFile *fRL = new TFile("out_RL.root");
	TFile *qq = TFile::Open("2f_Z.root");
	TFile *qqqq = TFile::Open("4f_WW.root");
	TFile *qqll = TFile::Open("4f_ZZ.root");
	TFile *qqlnu = TFile::Open("signal.root");
	
	TH1D* hqqtracks;
	TH1D* hqqtotalpt;
	TH1D* hqqtotalE;
	TH1D* hqqtotalM;
	TH1D* hqqyminus;
	TH1D* hqqyplus;

	TH1D* hqqqqtracks;
	TH1D* hqqqqtotalpt;
	TH1D* hqqqqtotalE;
	TH1D* hqqqqtotalM;
	TH1D* hqqqqyminus;
	TH1D* hqqqqyplus;

	TH1D* hqqlltracks;
	TH1D* hqqlltotalpt;
	TH1D* hqqlltotalE;
	TH1D* hqqlltotalM;
	TH1D* hqqllyminus;
	TH1D* hqqllyplus;

	TH1D* hstracks;
	TH1D* hstotalpt;
	TH1D* hstotalE;
	TH1D* hstotalM;
	TH1D* hsyminus;
	TH1D* hsyplus;
	

	

	qq->GetObject("htotaltracks", hqqtracks);
	qq->GetObject("htotalPt", hqqtotalpt);
	qq->GetObject("htotalE",hqqtotalE);
	qq->GetObject("htotalM",hqqtotalM);
	qq->GetObject("htotalym",hqqyminus);
	qq->GetObject("hyp",hqqyplus);
	
	qqqq->GetObject("htotaltracks", hqqqqtracks);
	qqqq->GetObject("htotalPt", hqqqqtotalpt);
	qqqq->GetObject("htotalE",hqqqqtotalE);
	qqqq->GetObject("htotalM",hqqqqtotalM);
	qqqq->GetObject("htotalym",hqqqqyminus);
	qqqq->GetObject("hyp",hqqqqyplus);

	qqll->GetObject("htotaltracks", hqqlltracks);
	qqll->GetObject("htotalPt", hqqlltotalpt);
	qqll->GetObject("htotalE",hqqlltotalE);
	qqll->GetObject("htotalM",hqqlltotalM);
	qqll->GetObject("htotalym",hqqllyminus);
	qqll->GetObject("hyp",hqqllyplus);

	qqlnu->GetObject("htotaltracks", hstracks);
	qqlnu->GetObject("htotalPt", hstotalpt);
	qqlnu->GetObject("htotalE",hstotalE);
	qqlnu->GetObject("htotalM",hstotalM);
	qqlnu->GetObject("htotalym",hsyminus);
	qqlnu->GetObject("hyp",hsyplus);
	
	hqqtracks->SetLineColor(kRed);
	hqqqqtracks->SetLineColor(kBlue);
	hqqlltracks->SetLineColor(kGreen);
	hstracks->SetLineColor(kBlack);

	hqqtotalpt->SetLineColor(kRed);
	hqqqqtotalpt->SetLineColor(kBlue);
	hqqlltotalpt->SetLineColor(kGreen);
	hstotalpt->SetLineColor(kBlack);

	hqqtotalE->SetLineColor(kRed);
	hqqqqtotalE->SetLineColor(kBlue);
	hqqlltotalE->SetLineColor(kGreen);
	hstotalE->SetLineColor(kBlack);

	hqqtotalM->SetLineColor(kRed);
	hqqqqtotalM->SetLineColor(kBlue);
	hqqlltotalM->SetLineColor(kGreen);
	hstotalM->SetLineColor(kBlack);

	hqqyminus->SetLineColor(kRed);
	hqqqqyminus->SetLineColor(kBlue);
	hqqllyminus->SetLineColor(kGreen);
	hsyminus->SetLineColor(kBlack);

	hqqyplus->SetLineColor(kRed);
	hqqqqyplus->SetLineColor(kBlue);
	hqqllyplus->SetLineColor(kGreen);
	hsyplus->SetLineColor(kBlack);


c1->SetTicks(1,1);


//PLOT TRACKS HERE /////////////  /*
/*
	hqqtracks->SetTitle(" ");
	hqqtracks->SetStats(kFALSE);
	hqqtracks->Draw();
	hqqqqtracks->Draw("same");
	hqqlltracks->Draw("same");
	hstracks->Draw("same");
	const float dx=0.36;
	const float dy=0.26;
	//Hfloat xlmin = 0.12;
	//float ylmin = 0.56;
TLegend* leg = new TLegend(xlmin,ylmin,xlmin+dx,ylmin+dy);
//leg->SetHeader("#sqrt{s} = 500 GeV: e^{-}_{L} e^{+}_{R} #rightarrow W^{-} W^{+} #rightarrow q q l #nu","C");
leg->AddEntry(hqqtracks,"ee #rightarrow qq","l");
leg->AddEntry(hqqqqtracks,"WW #rightarrow qqqq","l");
leg->AddEntry(hqqlltracks,"ZZ #rightarrow qqll","l");
leg->AddEntry(hstracks,"WW #rightarrow qql #nu","l");

leg->SetBorderSize(1);                          // Include border
leg->Draw();

*/


///////////////////////////
//PLOT pt HERE /////////////
/*
	hqqtotalpt->SetTitle(" ");
	hqqtotalpt->SetStats(kFALSE);
	hqqtotalpt->Draw();
	hqqqqtotalpt->Draw("same");
	hqqlltotalpt->Draw("same");
	hstotalpt->Draw("same");
	const float dx=0.36;
	const float dy=0.26;
	//Hfloat xlmin = 0.12;
	//float ylmin = 0.56;
TLegend* leg = new TLegend(xlmin,ylmin,xlmin+dx,ylmin+dy);
//leg->SetHeader("#sqrt{s} = 500 GeV: e^{-}_{L} e^{+}_{R} #rightarrow W^{-} W^{+} #rightarrow q q l #nu","C");
leg->AddEntry(hqqtotalpt,"ee #rightarrow qq","l");
leg->AddEntry(hqqqqtotalpt,"WW #rightarrow qqqq","l");
leg->AddEntry(hqqlltotalpt,"ZZ #rightarrow qqll","l");
leg->AddEntry(hstotalpt,"WW #rightarrow qql #nu","l");

leg->SetBorderSize(1);                          // Include border
leg->Draw();
*/

//PLOT E HERE /////////////
/*
	hqqtotalE->SetTitle(" ");
	hqqtotalE->SetStats(kFALSE);
	hqqtotalE->Draw();
	hqqqqtotalE->Draw("same");
	hqqlltotalE->Draw("same");
	hstotalE->Draw("same");
	const float dx=0.36;
	const float dy=0.26;
	//Hfloat xlmin = 0.12;
	//float ylmin = 0.56;
TLegend* leg = new TLegend(xlmin,ylmin,xlmin+dx,ylmin+dy);
//leg->SetHeader("#sqrt{s} = 500 GeV: e^{-}_{L} e^{+}_{R} #rightarrow W^{-} W^{+} #rightarrow q q l #nu","C");
leg->AddEntry(hqqtotalE,"ee #rightarrow qq","l");
leg->AddEntry(hqqqqtotalE,"WW #rightarrow qqqq","l");
leg->AddEntry(hqqlltotalE,"ZZ #rightarrow qqll","l");
leg->AddEntry(hstotalE,"WW #rightarrow qql #nu","l");

leg->SetBorderSize(1);                          // Include border
leg->Draw();
	

*/
/*
//PLOT M HERE /////////////
	hqqtotalM->SetTitle(" ");
	hqqtotalM->SetStats(kFALSE);
	hqqtotalM->Draw();
	hqqqqtotalM->Draw("same");
	hqqlltotalM->Draw("same");
	hstotalM->Draw("same");
	const float dx=0.36;
	const float dy=0.26;
	//Hfloat xlmin = 0.12;
	//float ylmin = 0.56;
TLegend* leg = new TLegend(xlmin,ylmin,xlmin+dx,ylmin+dy);
//leg->SetHeader("#sqrt{s} = 500 GeV: e^{-}_{L} e^{+}_{R} #rightarrow W^{-} W^{+} #rightarrow q q l #nu","C");
leg->AddEntry(hqqtotalM,"ee #rightarrow qq","l");
leg->AddEntry(hqqqqtotalM,"WW #rightarrow qqqq","l");
leg->AddEntry(hqqlltotalM,"ZZ #rightarrow qqll","l");
leg->AddEntry(hstotalM,"WW #rightarrow qql #nu","l");

leg->SetBorderSize(1);                          // Include border
leg->Draw();
*/
/*
//PLOT y- HERE /////////////
	hqqyminus->SetTitle(" ");
	hqqyminus->SetStats(kFALSE);
	hqqyminus->Draw();
	hqqqqyminus->Draw("same");
	hqqllyminus->Draw("same");
	hsyminus->Draw("same");
	const float dx=0.36;
	const float dy=0.26;
	//Hfloat xlmin = 0.12;
	//float ylmin = 0.56;
TLegend* leg = new TLegend(xlmin,ylmin,xlmin+dx,ylmin+dy);
//leg->SetHeader("#sqrt{s} = 500 GeV: e^{-}_{L} e^{+}_{R} #rightarrow W^{-} W^{+} #rightarrow q q l #nu","C");
leg->AddEntry(hqqyminus,"ee #rightarrow qq","l");
leg->AddEntry(hqqqqyminus,"WW #rightarrow qqqq","l");
leg->AddEntry(hqqllyminus,"ZZ #rightarrow qqll","l");
leg->AddEntry(hsyminus,"WW #rightarrow qql #nu","l");

leg->SetBorderSize(1);                          // Include border
leg->Draw();
*/
//PLOT y+ HERE /////////////

	hqqyplus->SetTitle(" ");
	hqqyplus->SetStats(kFALSE);
	hqqyplus->Draw();
	hqqqqyplus->Draw("same");
	hqqllyplus->Draw("same");
	hsyplus->Draw("same");
	const float dx=0.36;
	const float dy=0.26;
	//Hfloat xlmin = 0.12;
	//float ylmin = 0.56;
TLegend* leg = new TLegend(xlmin,ylmin,xlmin+dx,ylmin+dy);
//leg->SetHeader("#sqrt{s} = 500 GeV: e^{-}_{L} e^{+}_{R} #rightarrow W^{-} W^{+} #rightarrow q q l #nu","C");
leg->AddEntry(hqqyplus,"ee #rightarrow qq","l");
leg->AddEntry(hqqqqyplus,"WW #rightarrow qqqq","l");
leg->AddEntry(hqqllyplus,"ZZ #rightarrow qqll","l");
leg->AddEntry(hsyplus,"WW #rightarrow qql #nu","l");

leg->SetBorderSize(1);                          // Include border
leg->Draw();

//plot 1 of 4 vars

//TH1D * hLRutype = (TH1D*)fLR->Get(hutype.c_str());
//TH1D * hLRdtype = (TH1D*)fLR->Get(hdtype.c_str());
//TH1D * hLRltype = (TH1D*)fLR->Get(hltype.c_str());
//TH1D * hLRvtype = (TH1D*)fLR->Get(hvtype.c_str());

//TH1D * hRLutype = (TH1D*)fRL->Get(hutype.c_str());
//TH1D * hRLdtype = (TH1D*)fRL->Get(hdtype.c_str());
//TH1D * hRLltype = (TH1D*)fRL->Get(hltype.c_str());
//TH1D * hRLvtype = (TH1D*)fRL->Get(hvtype.c_str());

//hLRutype->GetYaxis()->SetTitleOffset(1.4);
//hLRutype->SetMaximum(ymax);
//hLRutype->SetStats(kFALSE);


/*hLRutype->SetLineColor(kRed);
hLRutype->Draw("ehist");
hLRdtype->SetLineColor(kBlue);
hLRdtype->Draw("same");
hLRltype->SetLineColor(kGreen);
hLRltype->Draw("same");
hLRvtype->SetLineColor(kBlack);
hLRvtype->Draw("same");
*/
// Also need a legend ...
/*const float dx=0.36;
const float dy=0.26;
float xlmin = 0.12;
float ylmin = 0.56; */
/*TLegend* leg = new TLegend(xlmin,ylmin,xlmin+dx,ylmin+dy);
leg->SetHeader("#sqrt{s} = 500 GeV: e^{-}_{L} e^{+}_{R} #rightarrow W^{-} W^{+} #rightarrow q q l #nu","C");
leg->AddEntry(hLRutype,"u-like q","l");
leg->AddEntry(hLRdtype,"d-like q","l");
leg->AddEntry(hLRltype,"lepton","l");
leg->AddEntry(hLRvtype,"#nu","l");

leg->SetBorderSize(1);                          // Include border
leg->Draw();
*/
//c1->Print("plot1.pdf").c_str());

}
