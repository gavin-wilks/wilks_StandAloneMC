/*******************************************************************
* Class for Monte Carlo Minimization                               *
* Log: Class started 29-07-2016                                    *
* Author: Y. Ayyad and W. Mittig (NSCL)                            *
********************************************************************/

#ifndef ATMCQMINIMIZATION_H
#define ATMCQMINIMIZATION_H

//ATTPCROOT
#include "ATMinimization.hh"
#include "ATHit.hh"
#include "ATDigiPar.hh"
#include "ATTrack.hh"
#include "AtTpcMap.h"

// FairRoot classes
#include "FairRuntimeDb.h"
#include "FairRun.h"

//ROOT
#include "TRotation.h"
#include "TMatrixD.h"
#include "TArrayD.h"
#include "TVector3.h"
#include "TH2Poly.h"
#include "TGraph.h"
#include "TRandom.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"

//System
#include <iostream>  //wm
#include <fstream>   // file stream wm
#include <cstdlib>   //wm
#include <vector>    //wm
#include <algorithm>

using std::vector;   //wm
using namespace std; //wm


#define cRED "\033[1;31m"
#define cYELLOW "\033[1;33m"
#define cNORMAL "\033[0m"
#define cGREEN "\033[1;32m"

class ATMCQMinimization : public ATMinimization{

      public:
	       ATMCQMinimization();
        ~ATMCQMinimization();

        void AddELossFunc(std::function<Double_t(Double_t,std::vector<Double_t>&)>& func);
        void AddRtoEFunc(std::function<Double_t(Double_t,std::vector<Double_t>&)>& func);
        void AddELossPar(std::vector<Double_t> (&par)[10]);
        void AddRtoEPar(std::vector<Double_t> (&par)[10]);
        void AddParticle(std::vector<std::pair<Int_t,Int_t>>& ptcl);
        void SetEntTB(Int_t value);
        void SetZGeoVertex(Bool_t value);
        void SetEntZ0(Double_t val);
        void SetBackWardPropagation(Bool_t value);
        void SetGainCalibration(Double_t value);
        void SetLongDiffCoef(Double_t value);
        void SetTranDiffCoef(Double_t value);
        void SetStepParameters(Double_t (&par)[10]);
        void SetRangeChi2(Bool_t value);
        static Double_t GetEloss(Double_t c0,std::vector<Double_t>& par);
	static Double_t GetEnergyFromRange(Double_t range,std::vector<Double_t>& par);


        Bool_t Minimize(Double_t* parameter,ATEvent *event);
        Bool_t MinimizeOpt(Double_t* parameter,ATEvent *event);
        Bool_t MinimizeOptMap(Double_t* parameter,ATEvent *event, TH2Poly* hPadPlane);
        Bool_t MinimizeOptMapAmp(Double_t* parameter,ATEvent *event, TH2Poly* hPadPlane,const multiarray& PadCoord);


        template <typename T, typename R>
        bool  MinimizeGen(Double_t* parameter,T* event ,const std::function<std::vector<R>*()>& func,TH2Poly* hPadPlane,const multiarray& PadCoord);


        std::vector<Double_t> GetPosXMin();
        std::vector<Double_t> GetPosYMin();
        std::vector<Double_t> GetPosZMin();
        std::vector<Double_t> GetPosXExp();
        std::vector<Double_t> GetPosYExp();
        std::vector<Double_t> GetPosZExp();
        std::vector<Double_t> GetPosXInt();
        std::vector<Double_t> GetPosYInt();
        std::vector<Double_t> GetPosZInt();
        std::vector<Double_t> GetPosXBack();
        std::vector<Double_t> GetPosYBack();
        std::vector<Double_t> GetPosZBack();
        Int_t GetMinimization();
        std::vector<ATHit> GetTBHitArray(Int_t TB,std::vector<ATHit> *harray);
        std::vector<std::function<Double_t(Double_t,std::vector<Double_t>&)>> *GetELossFunctionArray();

        void ResetParameters();




      protected:

       void CalibrateGain(std::vector<ATHit>* hitArray);
       void GetEnergy(Double_t M,Double_t IZ,Double_t BRO,Double_t &E);
       void GetBro(Double_t M,Double_t IZ,Double_t &BRO,Double_t E);
       Double_t GetSimThetaAngle(TVector3* pos, TVector3* posforw);
       void BackwardExtrapolation();
       void SetMap(AtTpcMap* map);
       Double_t GetChi2Pos(Int_t index,Int_t _iterCorrNorm,Int_t _par,
       Double_t *_xTBCorr,Double_t *_yTBCorr,Double_t *_zTBCorr);
       template <typename T>
       Double_t GetChi2Range(T* event,std::vector<Double_t> &_xTBCorr,std::vector<Double_t> &_yTBCorr,std::vector<Double_t> &_zTBCorr,Double_t sigma, Int_t npoints);

       Int_t GetTBHit(Int_t TB,std::vector<ATHit> *harray);
       TVector3 TransformIniPos(Double_t x,Double_t y, Double_t z); //Transforms initial position from Pad plane to Lab frame
       TVector3 InvTransIniPos(Double_t x,Double_t y, Double_t z); //Transforms lab frame to pad plane

       void PrintParameters(Int_t index);

       // New MC Amplitude functions
          void MCvar( double* parameter, int & MCmode, int & iconvar, double & x0MC, double & y0MC, double & z0MC,
                      double & thetaMC, double & phiMC, double & romin, double & x0MCv, double & y0MCv, double & z0MCv, 
		      double & thetaMCv, double & phiMCv, double & rominv);

          void QMCsim(double* parameter, double* Qsim, double *zsimq, double & QMCtotal,
                      double x0MC, double y0MC, double z0MC, double phiMCv, double thetaMCv,
                      double rominv, double & e0sm, multiarray PadCoord, TH2Poly *padplane);

          void Chi2MC(double Qtrack[10000], double ztrackq[10000], double & Qtracktotal,
                      double Qsim[10000], double zsimq[10000], double & QMCtotal,
                      double & Chi2fit, double & sigmaq, double & sigmaz);

             std::vector<ATHit> fHitTBArray;
             std::vector<ATHit> *fHitArray;

             Double_t fThetaMin;
             Double_t fEnerMin;
             TVector3 fPosMin;
             Double_t fBrhoMin;
             Double_t fBMin;
             Double_t fPhiMin;
             Double_t fDensMin;
             Double_t fVertexEner;
             TVector3 fVertexPos;
             Double_t fDens;
             Double_t fPressure;
             Double_t fGain;
             Double_t fMaxRange;
             Double_t fCoefL;
             Double_t fCoefT;



             std::vector<Double_t> fPosXmin;
             std::vector<Double_t> fPosYmin;
             std::vector<Double_t> fPosZmin;
             std::vector<Int_t>    fPosTBmin;
             std::vector<Double_t> fPosXexp;
             std::vector<Double_t> fPosYexp;
             std::vector<Double_t> fPosZexp;
             std::vector<Double_t> fPosXinter;
             std::vector<Double_t> fPosYinter;
             std::vector<Double_t> fPosZinter;
             std::vector<Int_t>    fPosTBinter;

             std::vector<Double_t> fPosXBack;
             std::vector<Double_t> fPosYBack;
             std::vector<Double_t> fPosZBack;

             std::vector<Double_t> fPosTBexp;
             std::vector<Double_t> fQmin;

             Double_t fDriftVelocity;
             Int_t fTBTime;
             Double_t fBField;
             Double_t fTiltAng;
             Double_t fThetaPad;
             Double_t fThetaLorentz;
             Double_t fThetaRot;
             Double_t fZk;
             Int_t fEntTB; //Beam entrance Time Bucket

             Double_t fEntZ0; // Calculated position of cathode/window
             Double_t fBeam_range; //Calculated range of the beam particle from range;

             //TRotation* fPadtoDetRot;

             //!AtTpcMap *fAtMapPtr;
             //!TH2Poly* fPadPlane;

             std::vector<std::function<Double_t(Double_t,std::vector<Double_t>&)>> fEloss_func_array; //!
             std::vector<std::function<Double_t(Double_t,std::vector<Double_t>&)>> fRtoE_func_array; //!
             std::vector<Double_t> fELossPar_array[10];
             std::vector<Double_t> fRtoEPar_array[10];
             std::vector<std::pair<Int_t,Int_t>> fParticleAZ;


             Bool_t kDebug;
             Bool_t kVerbose;
             Bool_t kPosChi2; //Enable the use of Position Chi2
             Bool_t kIsZGeoVertex; //Uses the relative Z vertex determined with the calibration performed with the original TB taken from parameter list
             Bool_t kBackWardProp; //Enables backward extrapolation if vertex is missing (default kTRUE)
             Bool_t kRangeChi2;

             //Global variables
             Double_t nucleons;
             Double_t m;
             Double_t dzstep;
             Int_t    integrationsteps;
             Double_t restmass;
             Double_t esm;
             Double_t protons;
             Double_t B0;
             Double_t B;

             std::vector<std::vector<ATHit>> *hitTBMatrix;//!
             Double_t *fXTBCorr;//!
             Double_t *fYTBCorr;//!
             Double_t *fZTBCorr;//!
             Int_t fIterCorrNorm;//!!

             Double_t fStep_par[10];
             Int_t fChi2Points;//!!


           ClassDef(ATMCQMinimization, 1);

};

#endif
