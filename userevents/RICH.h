// $Id: RICH.h,v 1.2 2021/04/05 12:59:34 ybedfer Exp ybedfer $

// $Log: RICH.h,v $
// Revision 1.2  2021/04/05 12:59:34  ybedfer
// "RICH_ESTIMATOR" changed from 7 (Max. Likelihood) to 10 (Ring Fit).
//
// Revision 1.1  2012/11/05 20:28:13  ybedfer
// Initial revision
//
// Revision 1.1  2004/02/11 14:36:10  ybedfer
// Initial revision
//

//#define RICH_CORRECT_THETA 4
//#define SCALE_RICH_INDEX

#ifndef RICH_h
#define RICH_h
/*!
  \class RICH
  \brief RICH 
 
  RICH gathers a number of RICH related methods:
  - PID methods: very basic ones, based on the Cerenkov angle theta. (Various
  determinations of the angle are returned by PaTrack::RichInf. One of
  them is selected by the macro "RICH_ESTIMATOR" defined infra. pID is based
  upon cuts on delta(theta) (w.r.t. to theoretical value): 1 cut for pID proper,
  and 1 extra cut for pi (in case of KID, else most relevant particle) veto.
  Cuts are momemtum (momentum-threshold in fact) dependent. Hence positive ID
  always implies above threshold. But no check is made that particle considered
  be  actually w/in RICH acceptance.)
  - Various other methods for:
    - retrieving thresholds,
    - correcting chi2 for diff between real index and production index,
    - etc...
  The current (i.e. corresponding to current run) indices have to be supplied
  via SetNRICH and SetNRICHProd.

  \author Yann.Bedfer@cern.ch
*/

#include <vector>
#include "PaTrack.h"

#define RICH_ESTIMATOR 10

class IndexData {
public:
  int run;
  double index;
  int time, length, nSpills;
  double prodIndex;
};

class RICH {
public:
  //! Constructor: Returned error has to be asserted
  RICH(double nSigmas, double nVeto);

  //! Returns pointer to RICH
  static RICH *Ptr() { return address; }

  // Refractive indices
  //! Accessor to index. Returns 0 by default.
  void GetNRICH(double &index, double &indAPV);
  //! Accessor to index used at production time.
  void GetNProd(double &index, double &indAPV);

  //! Retrieve indices (production one from setup and up-to-date one from DB), update RICH class member indices (and thresholds as well) and return them
  void UpdateIndices(int run, bool isMC,
		     double &runIndex,  double &prodIndex,
		     double &runIndAPV, double &prodIndAPV);
  void UpdateThresholds();

  //! Reset index
  void SetNRICH(double index, double indAPV);
  //! Reset prod. index
  void SetNRICHProd(double index, double indAPV);


  // PID
  //! eID: Return LS bit if eID, NLS bit if piVeto and NNLS bit if crude eID
  unsigned int eID(const PaTrack &trk);
  //! muID: Return LS bit if muID, NLS bit if piVeto and NNLS bit if crude muID
  unsigned int muID(const PaTrack &trk);
  //! KID: Return LS bit if KID and NLS bit if piVeto (N.B. other bits also set)
  unsigned int KID(const PaTrack &trk);
  //! piID: Return LS bit if piID and NLS bit if KVeto (N.B. other bits also set)
  unsigned int piID(const PaTrack &trk);
  //! pID: Return LS bit if pID and NLS bit if KVeto (N.B. other bits also set)
  unsigned int pID(const PaTrack &trk);

  //! (pi,K) ID's expressed as diffs (measure-expectation)/resol
  unsigned int piKID(const PaTrack &trk, double &piDth, double &KDth);
  //! Thresholds: are recalculated each time "CurrentIndex" changes.
  double eThr,    muThr,    piThr,    KThr,    pThr;
  double eThrAPV, muThrAPV, piThrAPV, KThrAPV, pThrAPV;

  // RICH studies
  //! Detailed pIDs for RICH studies.
  unsigned int eID(const PaTrack &trk, double &thRich, double &thMom);
  unsigned int KID(const PaTrack &trk, double &thRich, double &thMom);
  unsigned int piID(const PaTrack &trk, double &thRich, double &thMom);
  unsigned int pID(const PaTrack &trk, double &thRich, double &thMom);

  //! Return index of estimator used for theta Cerenkov in RichInf array
  int GetEstimator() { return RICH_ESTIMATOR; }

  //! Return Chi2's corrected for error on index
  bool CorrChiSq(double& chisqPi, double& chisqKa, double& chisqPr,
		 double  chisqRi, double theRi, double mom);

  //   ********** CURRENT INDEX VALUES **********
  //  There are 2 copies. Presently the one corresponding to the APV part is
  // not used for the thresholds determination. I am waiting for a method
  // allowing to determine whether a given track send its photons to the PMT or
  // APV part, to start  making use of it. In the mean time, the standard index
  // is always (for that, threshold, purpose) used. Note that it corresponds to
  // the UV region of the spectrum  in the 2002:2004 case and, in the case 2006
  // on, to the visible region, i.e. the most pervasive one in that case.
  //  On the other hand the Cerenkov angles returned by the RICH software are
  // renormalised so as to match the APV index => Therefore this APV index is
  // used for all other (than threshold) purposes.
  double CurrentIndex, CurrentIndAPV;
  double CurrentProdIndex, CurrentProdIndAPV;

private:
  static RICH *address;        //! (-) address of this object

  // pID
  //! Cut = dtheta(p)*nsigma; dtheta(p) = [0]+(x<[1])*[2]*(x-[1])^2; x = p-thr.
  double dthvsP[4], dthCut[4], dthVeto[4];

  //double M_pi, M2_pi, M_p, M2_p, M_K, M2_K;


};
// End of RICH class
#endif
