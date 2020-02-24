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

#include "ATEvent.hh"
#include "ATPad.hh"
#include "ATHit.hh"
#include "AtTpcMap.h"

#include "FairRootManager.h"
#include "FairLogger.h"
#include "FairRun.h"
#include "FairRunAna.h"

int target_thread_num = 4;


Int_t main()
{

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

    //Let's shoot in the dark, shall we?
    //Parameter Guesses
    parameter[0] = 0.;
    parameter[1] = 0.;
    parameter[2] = 0.;
    parameter[3] = 0.;
    parameter[4] = 0.;
    parameter[5] = 1.;
    parameter[6] = 4.0;
    parameter[7] = 1.;
   
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

  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  std::cout << std::endl << std::endl;
  std::cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << std::endl;
  std::cout << std::endl;
  return 0;

}
