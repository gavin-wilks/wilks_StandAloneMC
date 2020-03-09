#ifdef _OPENMP
#include <omp.h>
#endif

#include "MCSrc.hh"
#include "MCMinimization.hh"
#include "ATMCQMinimization.hh"
#include <ios>
#include <iostream>
#include <istream>
#include <limits>
#include <map>
#include <vector>

#include "TClonesArray.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreePlayer.h"
#include "TTreeReaderValue.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TH2Poly.h"
#include "TGraph.h"
#include "TApplication.h"

#include "ATEvent.hh"
#include "ATPad.hh"
#include "ATHit.hh"
#include "AtTpcMap.h"

#include "FairRootManager.h"
#include "FairLogger.h"
#include "FairRun.h"
#include "FairRunAna.h"

int target_thread_num = 4;


Int_t main(int argc, char* argv[])
{
    TApplication app("app",&argc,argv);
    gSystem->Load("libATTPCReco.so");
    //omp_set_num_threads(target_thread_num);

    TStopwatch timer;
    timer.Start();

    FairRunAna* run = new FairRunAna(); //Forcing a dummy run

    ATMCQMinimization *min = new ATMCQMinimization();
    min->ResetParameters();

    TString workdir = getenv("VMCWORKDIR");
    TString FileNameHead = "output";
    TString FilePath = workdir + "/macro/Unpack_HDF5/";
    TString FileNameTail = ".root";
    TString FileName     = FilePath + FileNameHead + FileNameTail;

    std::cout<<" Opening File : "<<FileName.Data()<<std::endl;
    TFile* file = new TFile(FileName.Data(),"READ");

    TTree* tree = (TTree*) file -> Get("cbmsim");
    Int_t nEvents = tree -> GetEntries();
    std::cout<<" Number of events : "<<nEvents<<std::endl;

    TTreeReader Reader1("cbmsim", file);
    TTreeReaderValue<TClonesArray> eventArray(Reader1, "ATEventH");
//    TTreeReaderValue<TClonesArray> houghArray(Reader1, "ATHough");

    Double_t* parameter = new Double_t[8];

    //Let's shoot in the dark
    //Parameter Guesses
    parameter[0] = 0.; //x
    parameter[1] = 0.; //y
    parameter[2] = 0.; //x
    parameter[3] = 0.; //TB
    parameter[4] = 0.; //phi
    parameter[5] = 1.; //romin
    parameter[6] = 0.; //theta
    parameter[7] = 1.; //hitID
   
    std::vector<Double_t> par[10];  //de/dx - E
    std::vector<Double_t> parRtoE[10]; // E - R
    Double_t stepPar[10];
    par[0] = {8.56,0.83,2.5,1.6,1.5,0.15,-1.0,-0.2,-0.17,-8.0,-0.4};
    parRtoE[0] = {0.63,-1.6,-1.0,0.5,19.0,-10.0,40.0};
    auto init = std::initializer_list<Double_t>({4.0,4.0,0.1,0.5,0.5,0.5,1,1,0,0.0});
    std::copy(init.begin(), init.end(),stepPar);
    std::vector<std::pair<Int_t,Int_t>> particle;
    particle.push_back(std::make_pair(4,2));

    min->AddELossPar(par);
    min->AddRtoEPar(parRtoE);
    min->SetStepParameters(stepPar);
    min->AddParticle(particle);
     
    std::function<Double_t(Double_t,std::vector<Double_t>&)> ELossFunc = std::bind(ATMCQMinimization::GetEloss,std::placeholders::_1,std::placeholders::_2);
    min->AddELossFunc(ELossFunc);

    std::function<Double_t(Double_t,std::vector<Double_t>&)> RtoEFunc = std::bind(ATMCQMinimization::GetEnergyFromRange,std::placeholders::_1,std::placeholders::_2);
    min->AddRtoEFunc(RtoEFunc);
    
    auto fAtMapPtr = new AtTpcMap();
    fAtMapPtr->GenerateATTPC();
    auto fPadPlane = fAtMapPtr->GetATTPCPlane(); 
    auto fAtPadCoord = fAtMapPtr->GetPadCoordArr();

    TCanvas * c1 = new TCanvas("c1", "c1",700,700);
    c1->Divide(2,2);
    TGraph * xy = new TGraph();
    TGraph * xz = new TGraph();
    TGraph * yz = new TGraph();

    nEvents = 100; //just show the first 100 events
      //#pragma omp parallel for ordered schedule(dynamic,1)
      for(Int_t i=0;i<nEvents;i++){
	  Reader1.Next();
          //while (Reader1.Next()) {
	      if(i == 15){
   	      std::cout<< "Event Number: " << i << std::endl;
              //Reader1.Next();

              ATEvent* event = (ATEvent*) eventArray->At(0);
              Int_t nHits = event->GetNumHits();
              std::vector<ATHit>* hitArray = event->GetHitArray();
              event->GetHitArrayObj();
              std::cout<<event->GetEventID()<<std::endl;
              hitArray->size();
             
  	      std::cout << "Hits: " << hitArray->size() << std::endl;

	      for(Int_t N=0;N<hitArray->size();N++){
              //std::cout << "I am finally working" << std::endl;
              ATHit hit = hitArray->at(N);
              TVector3 pos = hit.GetPosition();
              xy->SetPoint(N,pos.X(),pos.Y());
	      xz->SetPoint(N,pos.X(),pos.Z());
	      yz->SetPoint(N,pos.Y(),pos.Z());
             }

              /*std::vector<ATHit*>* hitbuff = new std::vector<ATHit*>;


                    for(Int_t iHit=0; iHit<nHits; iHit++){
                      ATHit hit = event->GetHit(iHit);
                      TVector3 hitPos = hit.GetPosition();
                      hitbuff->push_back(&hit);


                    }

		*/
              //std::cout<<hitbuff->size()<<std::endl;
              min->MinimizeOptMapAmp(&*parameter,&*event,&*fPadPlane,fAtPadCoord);
       	  }//if          
	}//loop

  //#pragma omp parallel for ordered schedule(dynamic,1)
  //for(Int_t i=0;i<100000;i++)std::cout<<" Hello ATTPCer! "<<std::endl;
  
  c1->cd(1);
  xy->SetMarkerSize(1);
  xy->SetMarkerStyle(4);
  xy->SetTitle("Y vs. X");
  xy->GetXaxis()->SetTitle("X");   
  xy->GetYaxis()->SetTitle("Y");
  xy->Draw("AP");

  c1->cd(2);
  xz->SetMarkerSize(1);
  xz->SetMarkerStyle(4);
  xz->SetTitle("Z vs. X");
  xz->GetXaxis()->SetTitle("X");   
  xz->GetYaxis()->SetTitle("Z");
  xz->Draw("AP");

  c1->cd(3);
  yz->SetMarkerSize(1);
  yz->SetMarkerStyle(4);
  yz->SetTitle("Z vs. Y");
  yz->GetXaxis()->SetTitle("Y");   
  yz->GetYaxis()->SetTitle("Z");
  yz->Draw("AP");
  
//c1->Draw();

  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  std::cout << std::endl << std::endl;
  std::cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << std::endl;
  std::cout << std::endl;

  app.Run();

  return 0;

}

