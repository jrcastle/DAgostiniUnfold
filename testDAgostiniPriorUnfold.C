#include "/home/j550c590/tdrstyle.C"
#include "DAgostiniUnfold.h"
#include "/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/interface/differentailBinning.h"
#include "/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfoldResponse.h"
#include "/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfoldBayes.h"
#include "/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfold.h"
#include "TLine.h"
#include "TH1D.h"
#include "TH2D.h" 
#include "TLegend.h"
#include "TLatex.h"
#include "TFile.h"
#include "TCanvas.h"

void testDAgostiniPriorUnfold(){

  bool validate = 1;
  bool dosys    = 0;
  int ptBin     = 0;
  int centBin   = 1;
  int VN        = 2;

  static const int NITER  = 4;
  static const int iter[] = {1, 2, 4, 8, 16, 32, 64, 128};
  int col[]               = {kOrange-2, kGreen+3, kCyan, kMagenta, kViolet-1, kBlue, kRed, kGray+2};


  setTDRStyle();
  TLatex latex;
  latex.SetNDC();

  TFile * fMC   = new TFile( Form("/rfs/jcastle/PbPb2015/mimicDataMC/cent%i-%i/pt%.2f-%.2f/MCTree.root", cent_min[centBin], cent_max[centBin], pt_min[ptBin], pt_max[ptBin]) );
  TFile * fData = new TFile( Form("/rfs/jcastle/PbPb2015/mimicDataMC/cent%i-%i/pt%.2f-%.2f/CastleEbyE.root", cent_min[centBin], cent_max[centBin], pt_min[ptBin], pt_max[ptBin]) );

  //-- CastleEbyE.root data members
  static const int NVn = 7;
  TH1D * hVnFull[NVn][NCENT_PBPB];
  TH2D * hResp[NVn][NCENT_PBPB];

  //-- Unfolding data
  TH2D * hresp[NCENT];
  TH1D * h1D[NCENT];
  TH1D * hTrue[NCENT];
  TH1D * hreco1[NCENT][NITER];
  TH1D * hreco1_RU[NCENT][NITER];
  TH1D * hratio_ME[NCENT][NITER];
  TH1D * hratio_RU[NCENT][NITER];
  RooUnfoldResponse * RUresp[NCENT];
  RooUnfoldBayes * RUBayes[NCENT][NITER];

  TLine * line = new TLine(0., 1., v2Max[ptBin][centBin], 1.);
  line->SetLineColor(2);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  TLegend * leg[NCENT][NITER];

  TCanvas * cIter[NCENT][NITER];

  //-- Get DATA
  for(int ivn = 1; ivn <=4; ivn++){
    if(ivn != VN) continue;
    for(int icent = 0; icent < NCENT_PBPB; icent++){
      if(icent == NCENT40_cutoff) break;
      hVnFull[ivn][icent]   = (TH1D*) fData->Get(Form("qwebye/hVnFull_%i_%i",ivn,icent));
      hResp[ivn][icent]     = (TH2D*) fData->Get(Form("qwebye/hResp_%i_%i",ivn,icent));
    }
  }

  //-- Condense from 2.5 to 10% cent binning
  for(int c = 0; c < NCENT; c++){

    if(c != centBin) continue;
    std::cout<<"!! Processing Cent = "<<c<<std::endl;

    //-- True distribution
    hTrue[c] = (TH1D*) fMC->Get("ebyeana/hV2True");
    hTrue[c]->SetLineColor(1);
    hTrue[c]->SetMarkerColor(1);
    hTrue[c]->SetMarkerStyle(21);

    //-- Observed distribution
    h1D[c] = (TH1D*)hVnFull[VN][4*c]->Clone(Form("h1D_%i_%i", VN, c));
    h1D[c]->Add(hVnFull[VN][4*c+1]);
    h1D[c]->Add(hVnFull[VN][4*c+2]);
    h1D[c]->Add(hVnFull[VN][4*c+3]);

    //-- Response matrix
    hresp[c] = (TH2D*) hResp[VN][4*c]->Clone(Form("hresp_%i_%i", VN, c));
    hresp[c]->SetOption("colz");
    hresp[c]->Add(hResp[VN][4*c+1]);
    hresp[c]->Add(hResp[VN][4*c+2]);
    hresp[c]->Add(hResp[VN][4*c+3]);

    RUresp[c] = new RooUnfoldResponse(0, 0, hresp[c], Form("response_%i_%i", VN, c));

    for(int i = 0; i<NITER; i++){

      if( !validate && i > 0 ) break;

      hreco1[c][i]    = 0;
      hreco1_RU[c][i] = 0;

      DAgostiniUnfold * unfolder = new DAgostiniUnfold(hresp[c]);
      if(dosys) unfolder->DoSystematics(); 
      unfolder->Unfold(h1D[c], iter[i], hTrue[c]);

      hreco1[c][i] = (TH1D*) unfolder->GetHReco( Form("hreco1_%ci_iter%i", c, iter[i]) );
      hreco1[c][i]->GetXaxis()->SetTitle("v_{2}");
      hreco1[c][i]->GetYaxis()->SetTitle("Events");
      hreco1[c][i]->SetLineColor(col[i]);
      hreco1[c][i]->SetMarkerColor(col[i]);
      hreco1[c][i]->SetMarkerStyle(20);

      delete unfolder;

      RUBayes[c][i] = new RooUnfoldBayes(RUresp[c], h1D[c], iter[i]);
      hreco1_RU[c][i] = (TH1D*) RUBayes[c][i]->Hreco();
      hreco1_RU[c][i]->GetXaxis()->SetTitle("v_{2}");
      hreco1_RU[c][i]->GetYaxis()->SetTitle("Events");
      hreco1_RU[c][i]->SetLineColor(col[i+4]);
      hreco1_RU[c][i]->SetMarkerColor(col[i+4]);
      hreco1_RU[c][i]->SetMarkerStyle(24);

      hratio_ME[c][i] = (TH1D*) hreco1[c][i]->Clone( Form("ratio_ME%i", i ) );
      hratio_ME[c][i]->Divide( hTrue[c] );
      hratio_ME[c][i]->GetXaxis()->SetTitle("v_{2}");
      hratio_ME[c][i]->GetYaxis()->SetTitle("Ratio");

      hratio_RU[c][i] = (TH1D*) hreco1_RU[c][i]->Clone( Form("ratio_RU%i", i ) );
      hratio_RU[c][i]->Divide( hTrue[c] );
      hratio_RU[c][i]->GetXaxis()->SetTitle("v_{2}");
      hratio_RU[c][i]->GetYaxis()->SetTitle("Ratio");

      leg[c][i] = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg[c][i]->SetFillStyle(0);
      leg[c][i]->SetBorderSize(0);
      leg[c][i]->AddEntry(hTrue[c], "MC Truth", "lp");
      leg[c][i]->AddEntry(hratio_ME[c][i], "CastleUnfold", "lp");
      leg[c][i]->AddEntry(hratio_RU[c][i], "RooUnfold", "lp");

      /*
      TH1D * hratioErr = (TH1D*) hratio->Clone("hratioErr");
      hratioErr->Reset();
      hratioErr->GetXaxis()->SetTitle("v_{2}");
      hratioErr->GetYaxis()->SetTitle("Error Bar Ratio");
      hratioErr->SetLineWidth(3);

      for(int i = 1; i <= hratio->GetNbinsX(); i++){
	
	double myErr = hreco1[centBin]->GetBinError(i);
	double ruErr = hreco1_RU[centBin]->GetBinError(i);
	if(ruErr == 0) continue;

	double ratio = myErr/ruErr;
	hratioErr->SetBinContent(i, ratio);
      }
      */

      cIter[c][i] = new TCanvas( Form("citer_c%i_iter%i", c, iter[i]), Form("citer_c%i_iter%i", c, iter[i]), 500, 1000);
      cIter[c][i]->Divide(1,2);
      cIter[c][i]->cd(1);
      cIter[c][i]->cd(1)->SetLogy();
      hTrue[c]->Draw();
      hreco1[c][i]->Draw("same");
      hreco1_RU[c][i]->Draw("same");
      leg[c][i]->Draw("same");
      latex.DrawLatex(0.2, 0.22, Form("N_{iter} = %i", iter[i]) );
      cIter[c][i]->cd(2);
      hratio_ME[c][i]->Draw();
      hratio_RU[c][i]->Draw("same");

      cIter[c][i]->SaveAs( Form("plots/citerPrior_c%i_iter%i.png", c, iter[i]) );

    } //-- End iter loop

  } //-- End cent loop

  /*
  TCanvas * c = new TCanvas("c","c", 500, 1000);
  c->Divide(1,2);
  c->cd(1);
  c->cd(1)->SetLogy();
  hreco1[centBin]->Draw();
  hreco1_RU[centBin]->Draw("same");

  c->cd(2);
  hratio->Draw();
  line->Draw("same");

  TCanvas * cc = new TCanvas("cc","cc",500, 500);
  cc->cd();
  hratioErr->Draw();
  line->Draw("same");
  */

}
