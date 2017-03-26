//--------------------------------------------------------------------------
/*
#include "Riostream.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>
*/
#include "TCanvas.h"
#include "TComplex.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TList.h"
#include "TMinuit.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"

TH2F *hxy_data;
TH2F *hxy_wmc;
TH2F *hxy_mc;
const Int_t rmc = 6400000;
Int_t nmc, nftpr;
Int_t iflnfun = 0, ntbin;
// Int_t ifnorm=0;
Double_t ndata, wNorm0, wNorm, ndata_fit, nmc_fit;
Int_t nhx = 18, nhy = 36;
Double_t rerr = 4.;
static Float_t MCX[rmc], MCY[rmc];

//______________________________________________________________________________
Double_t func1(Double_t x, Double_t y, Double_t *par) {
  Double_t value = 1. + par[0] * y + par[1] * y * y + par[2] * x * x;
  return value;
}

//______________________________________________________________________________
Double_t func2(Double_t x, Double_t y, Double_t *par) {
  Double_t value = 1. + 2. * par[0] * y +
                   (par[0] * par[0] + par[1] * par[1]) * y * y + par[2] * x * x;
  return value;
}

//______________________________________________________________________________
Double_t func3(Double_t x, Double_t y, Double_t *par) {
  Double_t value =
      1. + par[0] * y + par[1] * y * y + par[2] * x * x + par[3] * x * x * y;
  return value;
}

//______________________________________________________________________________
Double_t func4(Double_t x, Double_t y, Double_t *par) {
  Double_t value = 1. + par[0] * y + par[1] * y * y + par[2] * x * x +
                   par[3] * x * x * y + par[4] * x * x * x * x;
  return value;
}

//______________________________________________________________________________
Double_t func5(Double_t x, Double_t y, Double_t *par) {
  // Double_t value = 1.+ 2.*par[0]*y + (par[0]*par[0]+par[1]*par[1])*y*y +
  // par[2]*x*x + par[3]*x*x*y;
  // Double_t value = 1.+ par[0]*y + par[1]*y*y + par[2]*x*x + par[3]*fabs(x);
  Double_t value = 1. + par[0] * y + par[1] * y * y + par[2] * x * x +
                   par[3] * x * x * x * x;
  return value;
}
//______________________________________________________________________________
Double_t func6(Double_t x, Double_t y, Double_t *par) {
  Double_t value = par[3] * (1. + par[0] * y + par[1] * y * y + par[2] * x * x);
  return value;
}

//--------------------------------------------------------------------------

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
  Int_t iev, ix, iy;
  Double_t wMC, wMCa, wMCi, xMC, yMC, xFit, yFit;
  Double_t wMCt;
  Double_t wdat, wsim, edat, dat_sim, wsim0, etot2, esim;
  Double_t chisq = 0., chisqi, m2pi0, m2pi0sq, wcus;
  Int_t if_mchi2 = 1;

  hxy_wmc->Reset();
  for (iev = 0; iev < nmc; iev++) {
    xMC = MCX[iev];
    yMC = MCY[iev];
    if (iflnfun == 0)
      wMCi = func1(xMC, yMC, par);
    else if (iflnfun == 2)
      wMCi = func3(xMC, yMC, par);
    else if (iflnfun == 3)
      wMCi = func4(xMC, yMC, par);
    else if (iflnfun == 4)
      wMCi = func5(xMC, yMC, par);
    else if (iflnfun == 5)
      wMCi = func6(xMC, yMC, par);
    else
      wMCi = func2(xMC, yMC, par);
    hxy_wmc->Fill(xMC, yMC, wMCi);
  }

  wMCa = 0.;
  for (ix = 1; ix <= nhx; ix++) {
    for (iy = 1; iy <= nhy; iy++) {
      wdat = hxy_data->GetBinContent(ix, iy);
      edat = hxy_data->GetBinError(ix, iy);
      wsim0 = hxy_mc->GetBinContent(ix, iy);
      esim = hxy_mc->GetBinError(ix, iy);
      if (wsim0 > 1. && wdat > 1.) {
        if (wdat / edat < rerr || wsim0 / esim < rerr) {
          continue;
        }
      } else
        continue;
      wsim = hxy_wmc->GetBinContent(ix, iy);
      wMCa += wsim;
    }
  }
  wNorm = ndata_fit / wMCa;
  // printf("FCN:  wNorm %f, wMCa %f\n",wNorm,wMCa);

  for (ix = 1; ix <= nhx; ix++) {
    for (iy = 1; iy <= nhy; iy++) {
      wdat = hxy_data->GetBinContent(ix, iy);
      edat = hxy_data->GetBinError(ix, iy);
      wsim0 = hxy_mc->GetBinContent(ix, iy);
      if (wsim0 > 1. && wdat > 1.) {
        if (wdat / edat < rerr) {
          continue;
        }
      } else
        continue;
      wsim = hxy_wmc->GetBinContent(ix, iy);
      if (iflnfun != 5)
        wsim *= wNorm;
      etot2 = edat * edat;
      if (if_mchi2 != 0) {
        esim = hxy_wmc->GetBinError(ix, iy);
        if (iflnfun != 5)
          esim *= wNorm;
        etot2 += esim * esim;
      }
      dat_sim = wdat - wsim;
      chisqi = dat_sim * dat_sim / etot2;
      chisq += chisqi;
    }
  }
  if (nftpr == 3)
    printf("%d %f %f %f %f\n", ntbin, chisq, par[0], par[1], par[2]);
  if (nftpr == 4)
    printf("%d %f %f %f %f %f\n", ntbin, chisq, par[0], par[1], par[2], par[3]);
  if (nftpr == 5)
    printf("%d %f %f %f %f %f %f\n", ntbin, chisq, par[0], par[1], par[2],
           par[3], par[4]);
  f = chisq;
}

//--------------------------------------------------------------------------

void Accor_DPl_EtaprXY_2014_fit(Double_t a = -0.07, Double_t b = -0.07,
                                Double_t d = -0.06, Double_t we = 0.002,
                                Int_t iffit = 1, Int_t iflin = 0,
                                Int_t ifile = 0, Int_t lmc = 0, Int_t kfit = 2,
                                Double_t wup = 1., Double_t rer = 4.) {
  char *FileName, *FileName2, *FileName3, *FileName4, *FileName5, *FileName6;
  FileName = "DalPl_etapr_eta2pi0_data_accor_v2_read.dat";
  if (ifile == 1)
    FileName = "DalPl_etapr_eta2pi0_data_accor_v3_read.dat";

  Int_t nhi = 11;
  TH2F *Hi2D[nhi], *Hi2Da[nhi], *Hi2Db[nhi];
  TH1D *Hi1D[nhi], *Hi1Da[nhi], *Hi1Db[nhi];

  iflnfun = iflin;
  Int_t nfitpr = 3;
  if (iflnfun == 2)
    nfitpr = 4;
  if (iflnfun == 3)
    nfitpr = 5;
  if (iflnfun == 4)
    nfitpr = 4;
  if (iflnfun == 5)
    nfitpr = 4;
  nftpr = nfitpr;

  Int_t i;
  TMinuit *miniXY = new TMinuit(nfitpr + 1);
  miniXY->SetFCN(fcn);

  Double_t arglist[10];
  Int_t ierflg = 0;
  // arglist[0] = 1.;
  // miniXY->mnexcm("SET ERR", arglist ,1,ierflg);

  Double_t c = 0., rdtmc = 0.0865;
  const Int_t nprmax = 10;
  Double_t vstart[nprmax], step[nprmax];
  Int_t ipar = 0;
  vstart[ipar++] = a;
  vstart[ipar++] = b;
  vstart[ipar++] = d;
  if (nfitpr > 4)
    vstart[ipar++] = c;
  if (nfitpr > 5)
    vstart[ipar++] = 0.;
  if (iflnfun == 5)
    vstart[ipar++] = rdtmc;

  Double_t par4fit[nprmax], epar4fit[nprmax];
  char parname[2];
  for (i = 0; i < nfitpr; i++) {
    step[i] = we;
    sprintf(parname, "p%d", i);
    miniXY->mnparm(i, parname, vstart[i], step[i], 0, 0, ierflg);
    par4fit[i] = vstart[i];
  }
  if (iflnfun == 2 || iflnfun == 4) {
    arglist[0] = 1;
    arglist[1] = 2;
    arglist[2] = 3;
    miniXY->mnexcm("FIX", arglist, 3, ierflg);
  }

  Double_t rdn = 0.77, rup;
  rup = rdn + 0.4;
  Double_t xc = 1.1, yc = -0.85, yr = 0.8, xr;
  Double_t ymax = 1.15, ymin = -1., xmin = 0., xmax = 1.33;
  Int_t nhix = 18, nhiy = 36;

  Int_t iedk, idk, ndk;

  gROOT->SetStyle("Plain");

  TList *hl = (TList *)gDirectory->GetList();
  hl->Delete();

  TTree *fTree, *fTree2, *fTree3;
  Int_t Nentries, Nentries2, Nentries3, npromp = 0, nrndm = 0;
  Int_t iev, ix, iy;
  Double_t wMC, wMCa, wMCi, wMCt, wMCt2, xMC, yMC, xFit, yFit, x, y;
  Double_t wdat, wsim, edat, dat_sim, wsim0, etot2, esim;

  char HistName[120], HistTitle[120], XTitle[80], YTitle[80];
  char HistName2[120], HistName3[120], HistName4[120], HistName5[120];

  gStyle->SetStatW(0.47);
  gStyle->SetStatH(0.26);
  gStyle->SetStatY(1.);
  gStyle->SetTitleH(0.065);
  gStyle->SetTitleW(0.47);
  gStyle->SetOptStat(1000000);
  // gStyle->SetOptStat(10);
  gStyle->SetOptFit(111);
  // gStyle->SetStatFormat("10.0f");
  // gStyle->SetFitFormat("6.3g");

  TGaxis::SetMaxDigits(3);

  char cline[100];
  TLatex Lat;
  Lat.SetTextAlign(11);
  Lat.SetTextSize(0.085);

  TCanvas *c1 = new TCanvas("c1", "pict", 1, 1, 750, 650);
  c1->Divide(2, 2, 0.01, 0.01);
  Int_t ipad = 0;

  for (Int_t ip = 0; ip <= 2; ip++) {
    if (ip == 0) {
      FILE *f_dpl;
      f_dpl = fopen(FileName, "read");
      fscanf(f_dpl, "%d%lf%lf%d%lf%lf", &nhix, &xmin, &xmax, &nhiy, &ymin,
             &ymax);
      printf("Dal.pl: %d %f %f %d %f %f\n", nhix, xmin, xmax, nhiy, ymin, ymax);
      nhx = nhix;
      nhy = nhiy;
      hxy_data =
          new TH2F("hxy_data", "hxy_data", nhx, xmin, xmax, nhy, ymin, ymax);
      sprintf(HistTitle, "#eta'#rightarrow#eta2#pi^{0}#rightarrow6#gamma");
      sprintf(XTitle, "X");
      sprintf(YTitle, "Y");
      sprintf(cline, "(a)");

      ndk = nhx * nhy;
      for (idk = 1; idk <= ndk; idk++) {
        fscanf(f_dpl, "%d%d%lf%lf%lf%lf%lf%lf", &ix, &iy, &x, &y, &xMC, &yMC,
               &wdat, &edat);
        // printf("%d %d %f %f %f %f %f %f\n",ix,iy,x,y,xMC,yMC,wdat,edat);
        hxy_data->SetBinContent(ix, iy, wdat);
        hxy_data->SetBinError(ix, iy, edat);
      }
      fclose(f_dpl);
      Hi2D[ip] = (TH2F *)hxy_data->Clone();
      Hi2D[ip]->SetName("tmp1");
    }
    if (ip == 1) {
      hxy_mc = new TH2F("hxy_mc", "hxy_mc", nhx, xmin, xmax, nhy, ymin, ymax);
      hxy_wmc =
          new TH2F("hxy_wmc", "hxy_wmc", nhx, xmin, xmax, nhy, ymin, ymax);
      hxy_wmc->Sumw2();
      sprintf(HistTitle, "MC: #eta'#rightarrow#eta2#pi^{0}#rightarrow6#gamma");
      sprintf(cline, "(b)");
    }
    if (ip == 2) {
      sprintf(HistTitle, "WMC: #eta'#rightarrow#eta2#pi^{0}#rightarrow6#gamma");
      sprintf(cline, "(c)");
    }
    // delete gDirectory->FindObject(HistName);
    if (ip == 1) {
      nmc = 0;
      ndk = 32;
      Double_t eg, dk_x, dk_y, dk_mpi0eta[2], dk_m2pi0, dk_ax;
      char dk_fl_name[100];
      FILE *f_dk;
      for (idk = 1; idk <= ndk; idk++) {
        if (lmc > 0 && nmc > lmc)
          continue;
        sprintf(dk_fl_name, "dk_etapr2pi0eta2g_%d.dat", idk);
        f_dk = fopen(dk_fl_name, "read");
        for (iev = 0; iev < 200000; iev++) {
          fscanf(f_dk, "%d%lf", &iedk, &eg);
          fscanf(f_dk, "%lf%lf%lf%lf%lf", &dk_x, &dk_y, &dk_m2pi0,
                 &dk_mpi0eta[0], &dk_mpi0eta[1]);
          xMC = fabs(dk_x);
          yMC = dk_y;
          if (iflnfun == 0)
            wMCi = func1(xMC, yMC, par4fit);
          else if (iflnfun == 2)
            wMCi = func3(xMC, yMC, par4fit);
          else if (iflnfun == 3)
            wMCi = func4(xMC, yMC, par4fit);
          else if (iflnfun == 4)
            wMCi = func5(xMC, yMC, par4fit);
          else
            wMCi = func2(xMC, yMC, par4fit);
          hxy_mc->Fill(xMC, yMC);
          hxy_wmc->Fill(xMC, yMC, wMCi);
          MCX[nmc] = xMC;
          MCY[nmc] = yMC;
          nmc++;
        }
        fclose(f_dk);
      }

      printf("nMC: %d\n", nmc);

      ntbin = 0;
      ndata_fit = 0.;
      nmc_fit = 0.;
      wMCa = 0.;
      rerr = rer;
      for (ix = 1; ix <= nhx; ix++) {
        for (iy = 1; iy <= nhy; iy++) {
          wdat = hxy_data->GetBinContent(ix, iy);
          edat = hxy_data->GetBinError(ix, iy);
          wsim0 = hxy_mc->GetBinContent(ix, iy);
          if (wsim0 > 1. && wdat > 1.) {
            if (wdat / edat < rerr) {
              continue;
            }
          } else
            continue;
          ntbin++;
          ndata_fit += wdat;
          nmc_fit += wsim0;
          wsim = hxy_wmc->GetBinContent(ix, iy);
          wMCa += wsim;
        }
      }

      wMC = ndata_fit / nmc_fit;
      wNorm0 = wMC;
      printf("nData_fit: %f, nMC_fit: %f, wMC: %f, wMCa: %f\n", ndata_fit,
             nmc_fit, wMC, wMCa);

      Hi2D[ip] = (TH2F *)hxy_mc->Clone();
      Hi2D[ip]->SetName("tmp2");
    }
    if (ip == 2) {
      if (iffit == 0) {
        Hi2D[ip] = (TH2F *)hxy_wmc->Clone();
        Hi2D[ip]->SetName("tmp3");
      } else {
        arglist[0] = 0;
        arglist[1] = wup;
        if (kfit == 0)
          miniXY->mnexcm("CALL FCN", arglist, 2, ierflg);
        else if (kfit == 1) {
          arglist[1] *= 0.1;
          miniXY->mnexcm("SIMPLEX", arglist, 2, ierflg);
          // miniXY->mnexcm("SIMPLEX", arglist ,0,ierflg);
        } else {
          arglist[1] = wup;
          // miniXY->mnexcm("SET STRATEGY", arglist ,1,ierflg);
          miniXY->mnexcm("MIGRAD", arglist, 2, ierflg);
          // miniXY->mnexcm("MINOS", arglist ,2,ierflg);
          // miniXY->mnexcm("HESSE", arglist ,2,ierflg);
        }
        for (i = 0; i < nfitpr; i++) {
          miniXY->GetParameter(i, par4fit[i], epar4fit[i]);
          printf("param: %d, %f, %f\n", i, par4fit[i], epar4fit[i]);
        }
        // Print results
        /*
        Double_t amin,edm,errdef;
        Int_t nvpar,nparx,icstat;
        gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
        //gMinuit->mnprin(3,amin);
        */
        Hi2D[ip] = (TH2F *)hxy_wmc->Clone();
        Hi2D[ip]->SetName("tmp3");
      }
    }
    c1->cd(++ipad);
    // gPad->SetTopMargin(0.15);
    gPad->SetTopMargin(0.05);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.17);
    // gPad->SetRightMargin(0.03);
    gPad->SetRightMargin(0.12);

    // Hi2D[ip]->SetStats(1);
    Hi2D[ip]->SetStats(0);
    Hi2D[ip]->SetTitle("");
    Hi2D[ip]->GetXaxis()->SetLabelSize(0.06);
    Hi2D[ip]->GetYaxis()->SetLabelSize(0.06);
    Hi2D[ip]->GetXaxis()->SetTitleSize(0.07);
    Hi2D[ip]->GetYaxis()->SetTitleSize(0.065);
    Hi2D[ip]->GetXaxis()->SetTitle(XTitle);
    Hi2D[ip]->GetYaxis()->SetTitle(YTitle);
    Hi2D[ip]->GetXaxis()->CenterTitle();
    Hi2D[ip]->GetXaxis()->SetLabelOffset(0.01);
    Hi2D[ip]->GetYaxis()->SetLabelOffset(0.015);
    Hi2D[ip]->GetXaxis()->SetTitleOffset(1.1);
    Hi2D[ip]->GetYaxis()->SetTitleOffset(1.0);
    Hi2D[ip]->SetMinimum(0.);
    Hi2D[ip]->GetXaxis()->SetNdivisions(6);
    Hi2D[ip]->GetYaxis()->SetNdivisions(6);
    Hi2D[ip]->Draw("colz");
    Lat.SetTextSize(0.075);
    Lat.DrawLatex(xc, yc, cline);
  }
}

int main() { Accor_DPl_EtaprXY_2014_fit(); }
