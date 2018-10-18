void myplot(string htype="hx",float ymax=800.0, float xlmin=0.35, float ylmin=0.58){

string poln="LR";

TCanvas *c1 = new TCanvas("c1","multipads",800,600);
TFile *fLR = new TFile( ("out_"+poln+".root").c_str());
//TFile *fRL = new TFile("out_RL.root");

string hutype=htype+"u";
string hdtype=htype+"d";
string hltype=htype+"l";
string hvtype=htype+"v";

TH1D * hLRutype = (TH1D*)fLR->Get(hutype.c_str());
TH1D * hLRdtype = (TH1D*)fLR->Get(hdtype.c_str());
TH1D * hLRltype = (TH1D*)fLR->Get(hltype.c_str());
TH1D * hLRvtype = (TH1D*)fLR->Get(hvtype.c_str());

//TH1D * hRLutype = (TH1D*)fRL->Get(hutype.c_str());
//TH1D * hRLdtype = (TH1D*)fRL->Get(hdtype.c_str());
//TH1D * hRLltype = (TH1D*)fRL->Get(hltype.c_str());
//TH1D * hRLvtype = (TH1D*)fRL->Get(hvtype.c_str());

hLRutype->GetYaxis()->SetTitleOffset(1.4);
hLRutype->SetMaximum(ymax);
hLRutype->SetStats(kFALSE);

c1->SetTicks(1,1);

hLRutype->SetLineColor(kRed);
hLRutype->Draw("ehist");
hLRdtype->SetLineColor(kBlue);
hLRdtype->Draw("same");
hLRltype->SetLineColor(kGreen);
hLRltype->Draw("same");
hLRvtype->SetLineColor(kBlack);
hLRvtype->Draw("same");

// Also need a legend ...
const float dx=0.36;
const float dy=0.26;
//float xlmin = 0.12;
//float ylmin = 0.56;
TLegend* leg = new TLegend(xlmin,ylmin,xlmin+dx,ylmin+dy);
leg->SetHeader("#sqrt{s} = 500 GeV: e^{-}_{L} e^{+}_{R} #rightarrow W^{-} W^{+} #rightarrow q q l #nu","C");
leg->AddEntry(hLRutype,"u-like q","l");
leg->AddEntry(hLRdtype,"d-like q","l");
leg->AddEntry(hLRltype,"lepton","l");
leg->AddEntry(hLRvtype,"#nu","l");

leg->SetBorderSize(1);                          // Include border
leg->Draw();

c1->Print((htype+poln+".pdf").c_str());

}
