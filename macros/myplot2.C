void myplot2(string hutype="hsf", float ymax=500.0, bool logscale=false, float xlmin=0.48, float ylmin=0.61){

TCanvas *c1 = new TCanvas("c1","multipads",800,600);
TFile *fLR = new TFile("out_LR.root");
TFile *fRL = new TFile("out_RL.root");

//string hutype="hsf";

TH1D * hLRutype = (TH1D*)fLR->Get(hutype.c_str());
TH1D * hRLutype = (TH1D*)fRL->Get(hutype.c_str());

if(logscale)c1->SetLogy(1);
hLRutype->GetYaxis()->SetTitleOffset(1.4);
hLRutype->SetMaximum(ymax);
hLRutype->SetStats(kFALSE);

c1->SetTicks(1,1);

hLRutype->SetLineColor(kRed);
hLRutype->Draw("ehist");
hRLutype->SetLineColor(kBlue);
hRLutype->Draw("same");

// Also need a legend ...
const float dx=0.36;
const float dy=0.22;
TLegend* leg = new TLegend(xlmin,ylmin,xlmin+dx,ylmin+dy);
leg->SetHeader("#sqrt{s} = 500 GeV: e^{-} e^{+} #rightarrow W^{-} W^{+} #rightarrow q q l #nu","C");
leg->AddEntry(hLRutype,"e-(L), e+(R)","l");
leg->AddEntry(hRLutype,"e-(R), e+(L)","l");

leg->SetBorderSize(1);                          // Include border
leg->Draw();

if(logscale){
   c1->Print((hutype+"_Plot_log.pdf").c_str());
}
else{
   c1->Print((hutype+"_Plot.pdf").c_str());
}

}
