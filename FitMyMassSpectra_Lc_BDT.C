#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TInterpreter.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TPaveText.h>

#include "AliAODRecoDecayHF.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliHFInvMassFitter.h"
#endif


// MACRO to perform fits to D meson invariant mass spectra
// and store raw yields and cut object into a root output file
//
enum {kD0toKpi, kDplusKpipi};
enum {kBoth, kParticleOnly, kAntiParticleOnly};
enum {kExpo=0, kLinear, kPol2, kNoBk, kPow, kPowEx};
enum {kGaus=0, kDoubleGaus, kReflTempl};


// Common variables: to be configured by the user
Bool_t cutsappliedondistr=kFALSE;//kTRUE;
const Int_t nPtBins = 16;
const Double_t ptlims[nPtBins + 1] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};

const Double_t minRangeFit[]={2.2,   2.19,  2.2,  2.19, 2.19,2.19};//per hm
const Double_t maxRangeFit[]={2.38, 2.38, 2.38, 2.38, 2.38,2.38};//per hm

//pPb
const Double_t sigmaRef[]={
 0.00619594} ;
 ///  0.00619594,
// 0.00603889,//MC
// 0.00705208,
// 0.00763911,
// 0.00829985,
// 0.0095088
  
    
    
   //paper
//0.0061,    
//0.0062,
//0.0071,
//0.0075,
//0.0084,
//0.0095




const Double_t rebin[]={2,1,2,2,2,2};


Float_t sigmaFIT[nPtBins];//= {0.0078,0.0083,0.0090,0.010,0.012,0.012,0.013,0.013,0.014,0.015,0.0181,0.02258,0.0234,0.0286};
Float_t massRangeForCounting[]={0.03,0.03,0.03, 0.03, 0.035,0.05, 0.05, 0.05,0.055,0.055, 0.06,  0.08,  0.08,0.08};
Int_t typeb[14]={kPol2,kPol2,kPol2,kPol2,kPol2,kPol2,kPol2,kPol2,kPol2,kPol2,kPol2,kPol2,kPol2,kPol2};//kExpo,kExpo,kExpo,kExpo,kExpo,kExpo,kExpo,kExpo,kExpo,kExpo,kExpo,kExpo,kExpo};
Int_t types=kGaus;
Int_t optPartAntiPart=kBoth;//kAntiParticleOnly;//Both;
Int_t factor4refl=0;

TH2F* hPtMass=0x0;

int FitMassSpectra(Int_t analysisType=kD0toKpi,

 // TString fileNameb="AnalysisResults_1436_child1.root",
 // TString fileNamec="AnalysisResults_1436_child2.root",
//  TString fileNamed="AnalysisResults_1436_child3.root",
//  TString fileNamee="AnalysisResults_1436_child4.root",
                    
                    
       //             TString fileNameb="AnalysisResults_16qcent_1404.root",
       //             TString fileNamec="AnalysisResults_16qfast_1405.root",
      //              TString fileNamed="AnalysisResults_16tcent_1406.root",
       //             TString fileNamee="AnalysisResults_16tfast_1407.root",
//  AnalysisResults_1452_1_Lc_0100new.root
//AnalysisResults_1452_2_Lc_0100new.root
////AnalysisResults_1452_3_Lc_0100new.root
//AnalysisResults_1452_4_Lc_0100new.root






TString fileNameb="",//pass2
TString fileNamec="",
TString fileNamed="",
TString fileNamee="",
                    
                    
TString fileNamef="",
                    
		//    const char *CutsType="pp05V0A",
	//   const char *CutsType="pp010V0A",
//      coutputLcCutspp010V0A
		//   const char *CutsType="pp10_60V0A",
	//    const char *CutsType="pp60_100V0A",
	 //   const char *CutsType="ppProd_Cuts_40",
	    const char *CutsType="pp0_100V0A_c",
	    //const char *CutsType="pp0_100V0A",
		    const char *Centr="",
                    Int_t useTempl=0, //useTempl=1 if including reflections
                    bool isTHnSparse=false,
		    bool isARtask=false, 

   //  TString name="010_FixPass2MC"
   //  TString name="1060_FixPass2MC"
        TString name="60100_FixPass2MC"
   //    TString name="0100_FixtoPass2MC"
)

{
    //
    //_0_01V0M
    
    gROOT->SetStyle("Plain");
    gStyle->SetOptTitle(0);
    
    if(useTempl==1) types=kReflTempl;
    
    
    /*
     TObjArray* listFiles=new TObjArray();
     if(fileNameb!="") listFiles->AddLast(new TObjString(fileNameb.Data()));
     if(fileNamec!="") listFiles->AddLast(new TObjString(fileNamec.Data()));
     if(fileNamed!="") listFiles->AddLast(new TObjString(fileNamed.Data()));
     if(fileNamee!="") listFiles->AddLast(new TObjString(fileNamee.Data()));
     if(fileNamef!="") listFiles->AddLast(new TObjString(fileNamef.Data()));
     if(listFiles->GetEntries()==0){
     printf("Missing file names in input\n");
     return;
     }
     
     */
    
    
    
    TString fileName="output_final_kashmir.root";
    TFile* file = TFile::Open(fileName, "READ");
    
    
    TH1F** hmass=new TH1F*[nPtBins];
    TH1F** hTemplRefl=new TH1F*[nPtBins];
    TH1F** hSignMC=new TH1F*[nPtBins];
    
    TString filenameT="";
    TString filenameS="";
    
    Float_t massD;
    Bool_t retCode=kFALSE;
    
    Double_t nev=0;
    
        //
    for(Int_t i=0;i<nPtBins;i++) {
        hmass[i]=0x0;
        hTemplRefl[i]=0x0;
        hSignMC[i]=0x0;
    }
    
    
    // Loop over all pt bins
  /*  for(Int_t i=0;i<nPtBins;i++) {
        // Build the histogram name
        TString histName = Form("inv_mass_%d", i);
        
        // Get the histogram from the file
        hmass[i] = (TH1F*)file->Get(histName);
        
        // Check if histogram was properly loaded
        if(!hmass[i]){
            cout << "Histogram " << histName << " not found in file" << endl;
            continue;
        }*/
        
        
hmass[0] = (TH1F*)file->Get("inv_mass_10");
hmass[1] = (TH1F*)file->Get("inv_mass_15");
hmass[2] = (TH1F*)file->Get("inv_mass_20");
hmass[3] = (TH1F*)file->Get("inv_mass_25");
hmass[4] = (TH1F*)file->Get("inv_mass_30");
hmass[5] = (TH1F*)file->Get("inv_mass_35");
hmass[6] = (TH1F*)file->Get("inv_mass_40");
hmass[7] = (TH1F*)file->Get("inv_mass_45");
hmass[8] = (TH1F*)file->Get("inv_mass_50");
hmass[9] = (TH1F*)file->Get("inv_mass_55");
hmass[10] = (TH1F*)file->Get("inv_mass_60");
hmass[11] = (TH1F*)file->Get("inv_mass_65");
hmass[12] = (TH1F*)file->Get("inv_mass_70");
hmass[13] = (TH1F*)file->Get("inv_mass_75");
hmass[14] = (TH1F*)file->Get("inv_mass_80");
hmass[15] = (TH1F*)file->Get("inv_mass_85");
hmass[16] = (TH1F*)file->Get("inv_mass_90");
hmass[17] = (TH1F*)file->Get("inv_mass_95");
        
  

    
    
    massD=2.288;//TDatabasePDG::Instance()->GetParticle(4122)->Mass();//2.288
    //  if(!retCode){
    //   printf("ERROR in reading input files\n");
    
    //}
    
    
        
        Printf("massD------->%f",massD);
        
        TH1D* hCntSig1=new TH1D("hCntSig1","hCntSig1",nPtBins,ptlims);
        TH1D* hCntSig2=new TH1D("hCntSig2","hCntSig2",nPtBins,ptlims);
        TH1D* hNDiffCntSig1=new TH1D("hNDiffCntSig1","hNDiffCntSig1",nPtBins,ptlims);
        TH1D* hNDiffCntSig2=new TH1D("hNDiffCntSig2","hNDiffCntSig2",nPtBins,ptlims);
        TH1D* hSignal=new TH1D("hSignal","hSignal",nPtBins,ptlims);
        TH1D* hEv=new TH1D("hEv","Number of events",1,-0.5,0.5);
        TH1D* hRelErrSig=new TH1D("hRelErrSig","hRelErrSig",nPtBins,ptlims);
        TH1D* hInvSignif=new TH1D("hInvSignif","hInvSignif",nPtBins,ptlims);
        TH1D* hBackground=new TH1D("hBackground","hBackground",nPtBins,ptlims);
        TH1D* hBackgroundNormSigma=new TH1D("hBackgroundNormSigma","hBackgroundNormSigma",nPtBins,ptlims);
        TH1D* hSignificance=new TH1D("hSignificance","hSignificance",nPtBins,ptlims);
        TH1D* hMass=new TH1D("hMass","hMass",nPtBins,ptlims);
        TH1D* hSigma=new TH1D("hSigma","hSigma",nPtBins,ptlims);
        
        
        Int_t nMassBins=hmass[0]->GetNbinsX();
        Double_t hmin=hmass[0]->GetBinLowEdge(3);
        Double_t hmax=hmass[0]->GetBinLowEdge(nMassBins-2)+hmass[0]->GetBinWidth(nMassBins-2);
        Float_t minBinSum1=hmass[0]->FindBin(massD-massRangeForCounting[1]);
        Float_t maxBinSum1=hmass[0]->FindBin(massD+massRangeForCounting[1]);
        Int_t iPad=0;
        
        TF1** funBckStore1=new TF1*[nPtBins];
        TF1** funBckStore2=new TF1*[nPtBins];
        TF1** funBckStore3=new TF1*[nPtBins];
        
        
        AliHFInvMassFitter **fitter = new AliHFInvMassFitter*[nPtBins];
        Double_t arrchisquare[nPtBins];
        
        
        TCanvas* c1= new TCanvas("c1","MassSpectra");
        
        Int_t nx = 4, ny = 5;
        c1->Divide(nx, ny);
        
        TCanvas *myCanvas[nPtBins];
        TPaveText* ptBin[nPtBins];
        
        Double_t sig,errsig,s,errs,b,errb;
    
        
        
        //residuals
        TCanvas *cRes=new TCanvas("cRes","Residuales");
        Int_t nx_res = 4, ny_res = 5;
        cRes->Divide(nx_res, ny_res);
        
        
        
    
    
    for(Int_t iBin=0; iBin<nPtBins; iBin++){
        
        minBinSum1=hmass[iBin]->FindBin(massD-massRangeForCounting[0]);
        maxBinSum1=hmass[iBin]->FindBin(massD+massRangeForCounting[0]);
        
        c1->cd(iBin+1);
        //gPad->SetTicks();
        
        cout<<iBin<<"error finding"<<endl;
        Int_t origNbins=hmass[iBin]->GetNbinsX();
        hmass[iBin]->GetXaxis()->SetTitle("Invariant Mass (K#pi) (GeV/c^{2})");
        hmass[iBin]->GetYaxis()->SetTitle(Form("Entries / %1.0f MeV/c^{2}", (1000*(hmass[iBin]->GetXaxis()->GetBinWidth(1))*rebin[iBin])));
        hmass[iBin]->Rebin(rebin[iBin]);
        fitter[iBin]=new AliHFInvMassFitter(hmass[iBin],hmin,hmax,typeb[0],types);
        fitter[iBin]->SetUseLikelihoodFit();
        
        fitter[iBin]->SetRangeFit(minRangeFit[0],maxRangeFit[0]);
        //    rebin[iBin]=origNbins/fitter[iBin]->GetBinN();
        
        Printf("REBIN = %d", rebin[0]);
        
       // fitter[iBin]->SetInitialGaussianMean(massD);
           fitter[iBin]->SetFixGaussianMean(massD);
           fitter[iBin]->SetFixGaussianSigma(sigmaRef[0]);
        //if(iBin!=1)
        // fitter[iBin]->SetFixGaussianSigma(sigmaRef[iBin]);
       // fitter[iBin]->SetInitialGaussianSigma(0.00619594);
        
        
        Bool_t out=fitter[iBin]->MassFitter(kFALSE);
        cout<<"Bool: "<<out<<endl;
        
        if(!out) {
            fitter[iBin]->GetHistoClone()->Draw();
            continue;
        }
        cout<<"after check bool"<<endl;
        
        
        Double_t mass=fitter[iBin]->GetMean();
        Double_t massErr=fitter[iBin]->GetMeanUncertainty();
        Double_t sigma=fitter[iBin]->GetSigma();
        sigmaFIT[iBin]=fitter[iBin]->GetSigma();
        Double_t sigmaErr=fitter[iBin]->GetSigmaUncertainty();
        
        Float_t halfRange=3*sigma;
        Float_t minBinSum=hmass[iBin]->FindBin(massD-halfRange);
        Float_t maxBinSum=hmass[iBin]->FindBin(massD+halfRange);
        Printf("minBinSum (5sigma)=%f\t minBinSum (100MeV)=%f", minBinSum, minBinSum1);
        Printf("maxBinSum (5sigma)=%f\t maxBinSum (100MeV)=%f", maxBinSum, maxBinSum1);
        
        arrchisquare[iBin]=fitter[iBin]->GetReducedChiSquare();
        TF1* fB1=fitter[iBin]->GetBackgroundFullRangeFunc();
        TF1* fB2=fitter[iBin]->GetBackgroundRecalcFunc();
        TF1* fM=fitter[iBin]->GetMassFunc();
        funBckStore1[iBin]=(TF1*)fB1->Clone(Form("BkgFunction1_ptBin%d",iBin));
        funBckStore2[iBin]=(TF1*)fB2->Clone(Form("BkgFunction2_ptBin%d",iBin));
        funBckStore3[iBin]=(TF1*)fM->Clone(Form("MassFunction_ptBin%d",iBin));
        
        fitter[iBin]->DrawHere(gPad);
        fitter[iBin]->Signal(3,s,errs);
        fitter[iBin]->Background(3,b,errb);
        fitter[iBin]->Significance(3,sig,errsig);
        Double_t ry=fitter[iBin]->GetRawYield();
        Double_t ery=fitter[iBin]->GetRawYieldError();
        
        sigmaFIT[iBin]=fitter[iBin]->GetSigma();
        
        
        ptBin[iBin] = new TPaveText(0.6,0.55,0.8,0.6,"NDC");
        ptBin[iBin]->SetFillStyle(0);
        ptBin[iBin]->SetBorderSize(0);
        ptBin[iBin]->AddText(0.,0.,Form("%1.1f<p_{T}<%1.1f",ptlims[iBin],ptlims[iBin+1]));
        ptBin[iBin]->SetTextFont(42);
        ptBin[iBin]->SetTextSize(0.08);
        
        
        gPad->cd();
        ptBin[iBin]->Draw();
     //   gPad->Update();
        
        cRes->cd(iBin+1);
      //  gPad->SetTickx(2);
     //   gPad->SetTicky(2);
        fitter[iBin]->DrawHere(gPad);
        //
        //   cRes->cd(iBin+1+nPtBins);
        //   gPad->SetTickx(2);
        //subtraction background
        TH1F *hPulls = (TH1F*)hmass[iBin]->Clone("hPulls");
        TH1F *hResidualTrend = (TH1F*)hmass[iBin]->Clone("hResidualTrend");
        TH1F *hPullsTrend = (TH1F*)hmass[iBin]->Clone("hPullsTrend");
        TH1F *hsigsub = fitter[iBin]->GetOverBackgroundResidualsAndPulls(hPulls, hResidualTrend, hPullsTrend, minRangeFit[iBin],minRangeFit[iBin]);
            hResidualTrend->Draw();
        //fine subctration
        //if(nPtBins>1) {cMass[iMult]->cd(iPt+1);}
        //else {cMass[iMult]->cd(1);}
        //  fitter[iBin]->DrawHere(gPad);
        //
        //    myCanvas[iBin] = new TCanvas(Form("myCanvas_%d",iBin),Form("Invariant mass pt bin %d",iBin));
        //   fitter[iBin]->DrawHere(gPad);
        //   gPad->SetTicks();
        //   gPad->cd();
        //  ptBin[iBin]->Draw();
        //  gPad->Update();
        
        Float_t cntSig1=0.;
        Float_t cntSig2=0.;
        Float_t cntErr=0.;
        for(Int_t iMB=minBinSum1; iMB<=maxBinSum1; iMB++){
            Float_t bkg1=fB1 ? fB1->Eval(hmass[iBin]->GetBinCenter(iMB))/rebin[iBin] : 0;
            Float_t bkg2=fB2 ? fB2->Eval(hmass[iBin]->GetBinCenter(iMB))/rebin[iBin] : 0;
            cntSig1+=(hmass[iBin]->GetBinContent(iMB)-bkg1);
            cntSig2+=(hmass[iBin]->GetBinContent(iMB)-bkg2);
            cntErr+=(hmass[iBin]->GetBinContent(iMB));
        }
        hCntSig1->SetBinContent(iBin+1,cntSig1);
        hCntSig1->SetBinError(iBin+1,TMath::Sqrt(cntErr));
        // hNDiffCntSig1->SetBinContent(iBin+1,(s-cntSig1)/s);
        // hNDiffCntSig1->SetBinError(iBin+1,TMath::Sqrt(cntErr)/s);
        hNDiffCntSig1->SetBinContent(iBin+1,(ry-cntSig1)/ry);
        hNDiffCntSig1->SetBinError(iBin+1,TMath::Sqrt(cntErr)/ry);
        hCntSig2->SetBinContent(iBin+1,cntSig2);
        // hNDiffCntSig2->SetBinContent(iBin+1,(s-cntSig2)/s);
        // hNDiffCntSig2->SetBinError(iBin+1,TMath::Sqrt(cntErr)/s);
        hNDiffCntSig2->SetBinContent(iBin+1,(ry-cntSig2)/ry);
        hNDiffCntSig2->SetBinError(iBin+1,TMath::Sqrt(cntErr)/ry);
        hCntSig2->SetBinError(iBin+1,TMath::Sqrt(cntErr));
        // hSignal->SetBinContent(iBin+1,s);
        // hSignal->SetBinError(iBin+1,errs);
        hSignal->SetBinContent(iBin+1,ry);
        hSignal->SetBinError(iBin+1,ery);
        //    hEv->SetBinContent(1,nev);
        
        hRelErrSig->SetBinContent(iBin+1,errs/s);
        hInvSignif->SetBinContent(iBin+1,1/sig);
        hInvSignif->SetBinError(iBin+1,errsig/(sig*sig));
        hBackground->SetBinContent(iBin+1,b); //consider sigma
        hBackground->SetBinError(iBin+1,errb);
        hBackgroundNormSigma->SetBinContent(iBin+1,b/(3*fitter[iBin]->GetSigma())*(3*0.012)); //consider sigma
        hBackgroundNormSigma->SetBinError(iBin+1,errb);
        hSignificance->SetBinContent(iBin+1,sig);
        hSignificance->SetBinError(iBin+1,errsig);
        hMass->SetBinContent(iBin+1,mass);
        hMass->SetBinError(iBin+1,massErr);
        hSigma->SetBinContent(iBin+1,sigma);
        hSigma->SetBinError(iBin+1,sigmaErr);
        
        cout<<iBin<<"*%$"<<endl;
     //   if(iBin==0)return 0;
        
    }
    
    Float_t min=0.,max=0.,interval=0.;
  
    hEv->SetBinContent(1,nev);
    
    return 0;}
    
    
   
    

