#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TVector3.h"

void MakePlotNew(string sname="LR")
{
// Make histograms from elements of the TTree. Here the focus is on kinematic variables associated with the 4 fermions.

   TCanvas *c1 = new TCanvas("c1","multipads",800,600);

   // Open the file
   cout << sname << endl;
   string srootfile = "file_"+sname+".root";
   cout << srootfile << endl;
   TFile *f = TFile::Open(srootfile.c_str());
   f->Print();

   TH1D* hMud = new TH1D("hMud", "M(u,d)", 80, 40.0, 120.0);
   TH1D* hMlv = new TH1D("hMlv", "M(l,v)", 80, 40.0, 120.0);
   TH1D* hMave = new TH1D("hMave", "[M(u,d)+M(l,v)]/2", 80, 40.0, 120.0);
   TH1D* hcosud = new TH1D("hcosud", "cos(u,d)", 50, -1.0, 1.0);
   TH1D* hcosul = new TH1D("hcosul", "cos(u,l)", 50, -1.0, 1.0);
   TH1D* hcosdl = new TH1D("hcosdl", "cos(d,l)", 50, -1.0, 1.0);
   TH1D* hcoslv = new TH1D("hcoslv", "cos(l,v)", 50, -1.0, 1.0);
   TH1D* hcosljet = new TH1D("hcosljet", ";max{cos(l,j opening angle)};Events per bin", 50, -1.0, 1.0);
   TH1D* hxu = new TH1D("hxu", "; x(f) = E(f)/Ebeam; Events per bin", 50, 0.0, 1.0);
   TH1D* hxd = new TH1D("hxd", "; x(f) = E(f)/Ebeam; Events per bin", 50, 0.0, 1.0);
   TH1D* hxl = new TH1D("hxl", "; x(f) = E(f)/Ebeam; Events per bin", 50, 0.0, 1.0);
   TH1D* hxv = new TH1D("hxv", "; x(f) = E(f)/Ebeam; Events per bin", 50, 0.0, 1.0);
   TH1D* hxT = new TH1D("hxT", "; xT(v) = pT(v)/Ebeam; Events per bin", 50, 0.0, 1.0);
   TH1D* hsu = new TH1D("hsu", ";Signed sin#theta", 80, -1.0, 1.0);
   TH1D* hsd = new TH1D("hsd", ";Signed sin#theta", 80, -1.0, 1.0);
   TH1D* hsl = new TH1D("hsl", ";Signed sin#theta", 80, -1.0, 1.0);
   TH1D* hsv = new TH1D("hsv", ";Signed sin#theta", 80, -1.0, 1.0);
   TH1D* hsf = new TH1D("hsf", ";sin#theta (most forward visible fermion);Events per bin", 80, 0.0, 1.0);
   TH1D* hmqcosthl = new TH1D("hmqcosthl", "; -q cos#theta (l); Events per bin", 100, -1.0, 1.0);

   TTree *tree = (TTree *)f->Get("tree");

   TLorentzVector *f0 = 0;   tree->SetBranchAddress("MCf0",&f0);
   TLorentzVector *f1 = 0;   tree->SetBranchAddress("MCf1",&f1);
   TLorentzVector *f2 = 0;   tree->SetBranchAddress("MCf2",&f2);
   TLorentzVector *f3 = 0;   tree->SetBranchAddress("MCf3",&f3);
   int id0 = 0;              tree->SetBranchAddress("MCf0_PDG",&id0);
   int id1 = 0;              tree->SetBranchAddress("MCf1_PDG",&id1);
   int id2 = 0;              tree->SetBranchAddress("MCf2_PDG",&id2);
   int id3 = 0;              tree->SetBranchAddress("MCf3_PDG",&id3);

   int n = min(11400,int(tree->GetEntries()));
   cout << "n = " << n << endl;

   double Eb = 250.0;
   double leptonCharge;
   double acosthf;

   for (unsigned int ievt=0; ievt<n; ++ievt) {
      tree->GetEntry(ievt);
 //     cout << id0 << " " << id1 << " " << id2 << " " << id3 << endl;
      if(id2 == 11 || id2 == 13 || id2 == 15){
         leptonCharge = -1.0; 
      }else {
         leptonCharge =  1.0;
      }
      double cosang01 = cos(f0->Angle(f1->Vect()));
      double cosang02 = cos(f0->Angle(f2->Vect()));
      double cosang03 = cos(f0->Angle(f3->Vect()));
      double cosang12 = cos(f1->Angle(f2->Vect()));
      double cosang13 = cos(f1->Angle(f3->Vect()));
      double cosang23 = cos(f2->Angle(f3->Vect()));
      double cosanglj = max(cosang02,cosang12);

      std::pair < double, string > p01 = make_pair(cosang01,"ud");
      std::pair < double, string > p02 = make_pair(cosang02,"ul");
      std::pair < double, string > p12 = make_pair(cosang12,"dl");
      std::vector<pair <double, string> > vect;
      vect.push_back(p01);
      vect.push_back(p02);
      vect.push_back(p12);
      sort(vect.begin(), vect.end());    //Sorting based on first element in ascending order

      vector< pair <double, string>>::iterator it;
      for (it = vect.begin(); it<vect.end(); it++) {
          cout << (*it).first << " " << (*it).second << " " ;
      }
      cout << endl;

      cout << "1st " << vect[0].first << " " << vect[0].second << endl;
      cout << "2nd " << vect[1].first << " " << vect[1].second << endl;
      cout << "3rd " << vect[2].first << " " << vect[2].second << endl;
      cout << "Order: " << vect[0].second << " " << vect[1].second << " " << vect[2].second << endl;
      cout << " " << endl;


//      int icase = 0;
//      if(cosang01 > cos

      double xu = f0->E()/Eb;
      double xd = f1->E()/Eb;
      double xl = f2->E()/Eb;
      double xv = f3->E()/Eb;
      double xT = f3->Pt()/Eb;
      double costhu = f0->Pz()/f0->P();
      double costhd = f1->Pz()/f1->P();
      double costhl = f2->Pz()/f2->P();
      double costhv = f3->Pz()/f3->P();

      acosthf = 0.0;
      if(abs(costhu) > acosthf) acosthf = abs(costhu);
      if(abs(costhd) > acosthf) acosthf = abs(costhd);
      if(abs(costhl) > acosthf) acosthf = abs(costhl);

      double mqcosthl = -leptonCharge*costhl;
      if(ievt<4){
      cout << "EVENT " << ievt << endl;
      cout << f0->Px() << " " << f0->Py() << " " << f0->Pz() << " " << f0->E() << " " << f0->M() << endl;
      cout << f1->Px() << " " << f1->Py() << " " << f1->Pz() << " " << f1->E() << " " << f1->M() << endl;
      cout << f2->Px() << " " << f2->Py() << " " << f2->Pz() << " " << f2->E() << " " << f2->M() << endl;
      cout << f3->Px() << " " << f3->Py() << " " << f3->Pz() << " " << f3->E() << " " << f3->M() << endl;
      }

      hcosud->Fill(cosang01);
      hcosul->Fill(cosang02);
      hcosdl->Fill(cosang12);
      hcoslv->Fill(cosang23);
      hcosljet->Fill(cosanglj);
      hmqcosthl->Fill(mqcosthl);

      hxu->Fill(xu);
      hxd->Fill(xd);
      hxl->Fill(xl);
      hxv->Fill(xv);
      hxT->Fill(xT);

      hsu->Fill(copysign(1.0,costhu)*sin(acos(costhu)));
      hsd->Fill(copysign(1.0,costhd)*sin(acos(costhd)));
      hsl->Fill(copysign(1.0,costhl)*sin(acos(costhl)));
      hsv->Fill(copysign(1.0,costhv)*sin(acos(costhv)));
      hsf->Fill(sin(acos(acosthf)));

      TLorentzVector thadronic = *f0 + *f1;
      TLorentzVector tleptonic = *f2 + *f3;
      double mhad = thadronic.M();
      double mlv  = tleptonic.M();
      double mave = 0.5*(mhad+mlv);
      double mdif = mhad-mlv;
      //cout << "Mhad = " << mhad << " Mlv = " << mlv << endl;

      hMud->Fill(mhad);
      hMlv->Fill(mlv);
      hMave->Fill(mave);

/*         cout << "Angles " << endl;
         cout << "cos(u-d) " << cosang01 << endl;
         cout << "cos(l-v) " << cosang23 << endl;
         cout << "max cos(l-j) " << cosanglj << endl;
         cout << " " << endl; */
   }
   hcosud->Draw();
   hcosul->Draw();
   hcosdl->Draw();
   hcoslv->Draw();
   hcosljet->Draw();

   srootfile = "out_"+sname+".root";
   cout << "Writing output ROOT file " << srootfile << endl;
   TFile myfile(srootfile.c_str(),"RECREATE");

   hMud->Write();
   hMlv->Write();
   hMave->Write();
   hcosud->Write();
   hcosul->Write();
   hcosdl->Write();
   hcoslv->Write();
   hcosljet->Write();
   hmqcosthl->Write();
   hxu->Write();
   hxd->Write();
   hxl->Write();
   hxv->Write();
   hxT->Write();
   hsu->Write();
   hsd->Write();
   hsl->Write();
   hsv->Write();
   hsf->Write();
}
