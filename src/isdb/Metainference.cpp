/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* 
 The metainference update method was developed by Thomas Loehr
*/
#include "bias/Bias.h"
#include "bias/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "core/Value.h"
#include "tools/File.h"
#include "tools/OpenMP.h"
#include "tools/Random.h"
#include <cmath>
#include <ctime>

using namespace std;

namespace PLMD{
namespace isdb{

//+PLUMEDOC BIAS METAINFERENCE
/*
Calculate the Metainference energy for a set of back calculated experimental data.

The data can be averaged by using multiple replicas and weighted for a bias if present.
The functional form of Metainference can be chosen among four variants selected
with NOISE=GAUSS,MGAUSS,OUTLIERS,MOUTLIERS which correspond to modelling the noise for
the arguments as a single gaussian common to all the data points, a gaussian per data
point, a single long-tailed gaussian common to all the data points or a log-tailed
 gaussian per data point.

As from Metainference theory there are two sigma values: SIGMA_MEAN represent the
error of calculating an average quantity using a finite set of replica and should
be set as small as possible following the guidelines for replica-averaged simulations
in the framework of the Maximum Entropy Principle. Alternatively this can be obtained
automatically using the internal sigma mean optimisation (OPTSIGMAMEAN=SEM or FULL). 
SIGMA_BIAS is an uncertainty parameter, sampled by a MC algorithm in the bounded interval 
defined by SIGMA_MIN and SIGMA_MAX. The initial value is set at SIGMA0. The MC move is a 
random displacement of maximum value equal to DSIGMA. If the number of data point is
too large and the acceptance rate drops it is possible to make the MC move over mutually
exclusive, random subset of size MC_CHUNKSIZE and run more than one move setting MC_STRIDE
in such a way that MC_CHUNKSIZE*MC_STRIDE will cover all the data points.

Calculated and experimental data can be compared but for a scaling factor and/or an offset
using SCALEDATA and/or ADDOFFSET, the sampling is obtained by a MC algorithm either using
a flat or a gaussian prior setting it with SCALE_PRIOR or OFFSET_PRIOR.

\par Examples

In the following example we calculate a set of \ref RDC, take the replica-average of
them and comparing them with a set of experimental values. RDCs are compared with
the experimental data but for a multiplication factor SCALE that is also sampled by
MC on-the-fly

\verbatim
RDC ...
LABEL=rdc
SCALE=0.0001
GYROM=-72.5388
ATOMS1=22,23
ATOMS2=25,27
ATOMS3=29,31
ATOMS4=33,34
... RDC

METAINFERENCE ...
ARG=rdc.*
NOISETYPE=MGAUSS
PARAMETERS=1.9190,2.9190,3.9190,4.9190
SCALEDATA SCALE0=1 SCALE_MIN=0.1 SCALE_MAX=3 DSCALE=0.01 
SIGMA0=0.01 SIGMA_MIN=0.00001 SIGMA_MAX=3 DSIGMA=0.01 
SIGMA_MEAN=0.001
LABEL=spe
... METAINFERENCE 

PRINT ARG=spe.bias FILE=BIAS STRIDE=1 
\endverbatim

in the following example instead of using one uncertainty parameter per data point we use
a single uncertainty value in a long-tailed gaussian to take into account for outliers, furthermore
the data are weighted for the bias applied to other variables of the system.

\verbatim
cv1: TORSION ATOMS=1,2,3,4
cv2: TORSION ATOMS=2,3,4,5
mm: METAD ARG=cv1,cv2 HEIGHT=0.5 SIGMA=0.3,0.3 PACE=200 BIASFACTOR=8 WALKERS_MPI

METAINFERENCE ...
ARG=rdc.*,mm.bias
REWEIGHT
NOISETYPE=OUTLIERS
PARAMETERS=1.9190,2.9190,3.9190,4.9190
SCALEDATA SCALE0=1 SCALE_MIN=0.1 SCALE_MAX=3 DSCALE=0.01 
SIGMA0=0.01 SIGMA_MIN=0.00001 SIGMA_MAX=3 DSIGMA=0.01 
SIGMA_MEAN=0.001
LABEL=spe
... METAINFERENCE 
\endverbatim

(See also \ref RDC, \ref PBMETAD).

*/
//+ENDPLUMEDOC

class Metainference : public bias::Bias
{
  // experimental values
  vector<double> parameters;
  // noise type
  unsigned noise_type_;
  enum { GAUSS, MGAUSS, OUTLIERS, MOUTLIERS, GENERIC };
  // scale is data scaling factor
  // noise type
  unsigned scale_prior_;
  enum { SC_GAUSS, SC_FLAT };
  bool   doscale_;
  double scale_;
  double scale_mu_;
  double scale_min_;
  double scale_max_;
  double Dscale_;
  // scale is data scaling factor
  // noise type
  unsigned offset_prior_;
  bool   dooffset_;
  double offset_;
  double offset_mu_;
  double offset_min_;
  double offset_max_;
  double Doffset_;
  // sigma is data uncertainty
  vector<double> sigma_;
  vector<double> sigma_min_;
  vector<double> sigma_max_;
  vector<double> Dsigma_;
  // sigma_mean is uncertainty in the mean estimate
  vector<double> sigma_mean2_;
  vector<double> ftilde_;
  double Dftilde_;

  // sigma_mean rescue params
  double sigma_mean_correction_;
  double sm_mod_;
  double sm_mod_min_;
  double sm_mod_max_;
  double Dsm_mod_;
  double max_force_;

  // temperature in kbt
  double   kbt_;

  // Monte Carlo stuff
  vector<Random> random;
  unsigned MCsteps_;
  unsigned MCstride_;
  long unsigned MCaccept_;
  long unsigned MCacceptScale_;
  long unsigned MCacceptFT_;
  long unsigned MCtrial_;
  unsigned MCchunksize_;

  // output
  Value*   valueScale;
  Value*   valueOffset;
  Value*   valueAccept;
  Value*   valueAcceptScale;
  vector<Value*> valueSigma;
  vector<Value*> valueSigmaMean;
  vector<Value*> valueFtilde;
  Value*   valueSMmod;
  Value*   valueMaxForceMD;

  // restart
  unsigned write_stride_;
  OFile    sfile_;

  // others
  bool     firsttime;
  bool     master;
  bool     do_reweight;
  unsigned do_optsigmamean_;
  unsigned nrep_;
  unsigned replica_;
  unsigned narg;

  // we need this for the forces
  Atoms& atoms;

  double getEnergyMGGEN(const vector<double> &mean, const vector<double> &ftilde, const vector<double> &sigma, 
                        const double scale, const double offset, const double modifier);
  double getEnergySP(const vector<double> &mean, const vector<double> &sigma, 
                     const double scale, const double offset, const double modifier);
  double getEnergySPE(const vector<double> &mean, const vector<double> &sigma, 
                      const double scale, const double offset, const double modifier);
  double getEnergyGJ(const vector<double> &mean, const vector<double> &sigma, 
                     const double scale, const double offset, const double modifier);
  double getEnergyGJE(const vector<double> &mean, const vector<double> &sigma, 
                      const double scale, const double offset, const double modifier);
  void   doMonteCarlo(const vector<double> &mean, const double modifier);
  double getEnergyForceSP(const vector<double> &mean, const vector<double> &dmean_x, const vector<double> &dmean_b, const double modifier);
  double getEnergyForceSPE(const vector<double> &mean, const vector<double> &dmean_x, const vector<double> &dmean_b, const double modifier);
  double getEnergyForceGJ(const vector<double> &mean, const vector<double> &dmean_x, const vector<double> &dmean_b, const double modifier);
  double getEnergyForceGJE(const vector<double> &mean, const vector<double> &dmean_x, const vector<double> &dmean_b, 
                           const vector<double> &dsigma_mean2_x, const vector<double> &dsigma_mean2_b, const double modifier);
  double getEnergyForceMGGEN(const vector<double> &mean, const vector<double> &dmean_x, const vector<double> &dmean_b, 
                             const vector<double> &dsigma_mean2_x, const vector<double> &dsigma_mean2_b, const double modifier);
  void   writeStatus();
  
public:
  explicit Metainference(const ActionOptions&);
  ~Metainference();
  void calculate();
  void update();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Metainference,"METAINFERENCE")

void Metainference::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("optional","PARARG","reference values for the experimental data, these can be provided as arguments without derivatives"); 
  keys.add("optional","PARAMETERS","reference values for the experimental data");
  keys.add("compulsory","NOISETYPE","functional form of the noise (GAUSS,MGAUSS,OUTLIERS,MOUTLIERS)");
  keys.addFlag("REWEIGHT",false,"simple REWEIGHT using the latest ARG as energy"); 
  keys.addFlag("SCALEDATA",false,"Set to TRUE if you want to sample a scaling factor common to all values and replicas");  
  keys.addFlag("ADDOFFSET",false,"Set to TRUE if you want to sample an offset common to all values and replicas");  
  keys.add("compulsory","SCALE0","initial value of the scaling factor");
  keys.add("compulsory","SCALE_PRIOR","FLAT","either FLAT or GAUSSIAN");
  keys.add("optional","SCALE_MIN","minimum value of the scaling factor");
  keys.add("optional","SCALE_MAX","maximum value of the scaling factor");
  keys.add("optional","DSCALE","maximum MC move of the scaling factor");
  keys.add("compulsory","OFFSET0","initial value of the offset");
  keys.add("compulsory","OFFSET_PRIOR","FLAT","either FLAT or GAUSSIAN");
  keys.add("optional","OFFSET_MIN","minimum value of the offset");
  keys.add("optional","OFFSET_MAX","maximum value of the offset");
  keys.add("optional","DOFFSET","maximum MC move of the offset");
  keys.add("compulsory","SIGMA0","initial value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MIN","minimum value of the uncertainty parameter");
  keys.add("compulsory","SIGMA_MAX","maximum value of the uncertainty parameter");
  keys.add("optional","DSIGMA","maximum MC move of the uncertainty parameter");
  keys.add("optional","DFTILDE","maximum MC move of the scaling factor");
  keys.add("compulsory","OPTSIGMAMEAN","NONE","Set to NONE/SEM/FULL to manually set sigma mean, or to estimate it on the fly and use the safety check on forces");  
  keys.add("optional","SIGMA_MEAN0","starting value for the uncertainty in the mean estimate");
  keys.add("optional","SIGMA_MEAN_MOD0","starting value for sm modifier");
  keys.add("optional","SIGMA_MEAN_MOD_MIN","starting value for sm modifier");
  keys.add("optional","SIGMA_MEAN_MOD_MAX","starting value for sm modifier");
  keys.add("optional","DSIGMA_MEAN_MOD","step value for sm modifier");
  keys.add("optional","SIGMA_MEAN_CORRECTION","sigma_mean correction modifier");
  keys.add("optional","MAX_FORCE","maximum allowable force");
  keys.add("optional","TEMP","the system temperature - this is only needed if code doesnt' pass the temperature to plumed");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  keys.add("optional","MC_CHUNKSIZE","MC chunksize");
  keys.add("optional","STATUS_FILE","write a file with all the data usefull for restart/continuation of Metainference");
  keys.add("compulsory","WRITE_STRIDE","write the status to a file every N steps, this can be used for restart/continuation");
  keys.use("RESTART");
  useCustomisableComponents(keys);
  keys.addOutputComponent("sigma",        "default",      "uncertainty parameter");
  keys.addOutputComponent("sigmaMean",    "default",      "uncertainty in the mean estimate");
  keys.addOutputComponent("acceptSigma",  "default",      "MC acceptance");
  keys.addOutputComponent("acceptScale",  "SCALEDATA",    "MC acceptance");
  keys.addOutputComponent("weight",       "REWEIGHT",     "weights of the weighted average");
  keys.addOutputComponent("MetaDf",       "REWEIGHT",     "force on metadynamics");
  keys.addOutputComponent("scale",        "SCALEDATA",    "scale parameter");
  keys.addOutputComponent("offset",       "ADDOFFSET",    "offset parameter");
  keys.addOutputComponent("maxForceMD",   "OPTSIGMAMEAN", "max force on atoms");
  keys.addOutputComponent("smMod",        "OPTSIGMAMEAN", "modifier for all sigma mean");
}

Metainference::Metainference(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao), 
doscale_(false),
scale_mu_(0),
scale_min_(1),
scale_max_(-1),
Dscale_(-1),
dooffset_(false),
offset_mu_(0),
offset_min_(1),
offset_max_(-1),
Doffset_(-1),
firsttime(true),
sigma_mean_correction_(1.),
sm_mod_(1.),
random(3),
MCsteps_(1), 
MCstride_(1), 
MCaccept_(0), 
MCacceptScale_(0), 
MCtrial_(0),
MCchunksize_(0),
write_stride_(0),
do_reweight(false),
do_optsigmamean_(0),
atoms(plumed.getAtoms())
{
  // set up replica stuff 
  master = (comm.Get_rank()==0);
  if(master) {
    nrep_    = multi_sim_comm.Get_size();
    replica_ = multi_sim_comm.Get_rank();
  } else {
    nrep_    = 0;
    replica_ = 0;
  }
  comm.Sum(&nrep_,1);
  comm.Sum(&replica_,1);

  // reweight implies a different number of arguments (the latest one must always be the bias) 
  parseFlag("REWEIGHT", do_reweight);
  if(do_reweight&&nrep_==0) error("REWEIGHT can only be used in parallel with 2 or more replicas"); 
  narg = getNumberOfArguments();
  if(do_reweight) narg--;

  parseVector("PARAMETERS",parameters);
  if(parameters.size()!=static_cast<unsigned>(narg)&&!parameters.empty())
    error("Size of PARAMETERS array should be either 0 or the same as of the number of arguments in ARG1");

  vector<Value*> arg2;
  parseArgumentList("PARARG",arg2);
  if(!arg2.empty()) {
    if(parameters.size()>0) error("It is not possible to use PARARG and PARAMETERS together");
    if(arg2.size()!=narg) error("Size of PARARG array should be the same as number for arguments in ARG");
    for(unsigned i=0;i<arg2.size();i++){
      parameters.push_back(arg2[i]->get()); 
      if(arg2[i]->hasDerivatives()==true) error("PARARG can only accept arguments without derivatives");
    }
  }

  if(parameters.size()!=narg) 
    error("PARARG or PARAMETERS arrays should include the same number of elements as the arguments in ARG");

  string stringa_noise;
  parse("NOISETYPE",stringa_noise);
  if(stringa_noise=="GAUSS")           noise_type_ = GAUSS; 
  else if(stringa_noise=="MGAUSS")     noise_type_ = MGAUSS;
  else if(stringa_noise=="OUTLIERS")   noise_type_ = OUTLIERS;
  else if(stringa_noise=="MOUTLIERS")  noise_type_ = MOUTLIERS;
  else if(stringa_noise=="GENERIC")    noise_type_ = GENERIC;
  else error("Unknown noise type!"); 
 
  parse("DFTILDE",Dftilde_);
  ftilde_.resize(narg,0.);

  parse("WRITE_STRIDE",write_stride_);
  string status_file_name_;
  parse("STATUS_FILE",status_file_name_);
  if(status_file_name_=="") status_file_name_ = "MISTATUS"+getLabel();

  string stringa_optsigma;
  parse("OPTSIGMAMEAN", stringa_optsigma);
  if(stringa_optsigma=="NONE")      do_optsigmamean_=0;
  else if(stringa_optsigma=="SEM")  do_optsigmamean_=1;
  else if(stringa_optsigma=="FULL") do_optsigmamean_=2;

  parseFlag("SCALEDATA", doscale_);
  if(doscale_) {
    string stringa_noise;
    parse("SCALE_PRIOR",stringa_noise);
    if(stringa_noise=="GAUSSIAN")  scale_prior_ = SC_GAUSS; 
    else if(stringa_noise=="FLAT") scale_prior_ = SC_FLAT;
    else error("Unknown SCALE_PRIOR type!");
    parse("SCALE0",scale_);
    parse("DSCALE",Dscale_);
    if(Dscale_<0.) error("DSCALE must be set when using SCALEDATA");
    if(scale_prior_==SC_GAUSS) {
      scale_mu_=scale_;
    } else {
      parse("SCALE_MIN",scale_min_);
      parse("SCALE_MAX",scale_max_);
      if(scale_max_<scale_min_) error("SCALE_MAX must be set and be larger than SCALE_MIN when using SCALE_PRIOR=FLAT");
    }
  } else {
    scale_=1.0;
  }

  parseFlag("ADDOFFSET", dooffset_);
  if(dooffset_) {
    string stringa_noise;
    parse("OFFSET_PRIOR",stringa_noise);
    if(stringa_noise=="GAUSSIAN")  offset_prior_ = SC_GAUSS; 
    else if(stringa_noise=="FLAT") offset_prior_ = SC_FLAT;
    else error("Unknown OFFSET_PRIOR type!");
    parse("OFFSET0",offset_);
    parse("DOFFSET",Doffset_);
    if(offset_prior_==SC_GAUSS) {
      offset_mu_=offset_;
      if(Doffset_<0.) error("DOFFSET must be set when using OFFSET_PRIOR=GAUSS");
    } else {
      parse("OFFSET_MIN",offset_min_);
      parse("OFFSET_MAX",offset_max_);
      if(Doffset_<0) Doffset_ = 0.05*(offset_max_ - offset_min_);
      if(offset_max_<offset_min_) error("OFFSET_MAX and OFFSET_MIN must be set when using OFFSET_PRIOR=FLAT");
    }
  } else {
    offset_=0.;
  }

  vector<double> readsigma;
  parseVector("SIGMA0",readsigma);
  if((noise_type_!=MGAUSS&&noise_type_!=MOUTLIERS&&noise_type_!=GENERIC)&&readsigma.size()>1) 
    error("If you want to use more than one SIGMA you should use NOISETYPE=MGAUSS|MOUTLIERS");
  if(noise_type_==MGAUSS||noise_type_==MOUTLIERS||noise_type_==GENERIC) {
    if(readsigma.size()==narg) {
      sigma_.resize(narg);
      sigma_=readsigma;
    } else if(readsigma.size()==1) {
      sigma_.resize(narg,readsigma[0]);
    } else {
      error("SIGMA0 can accept either one single value or as many values as the number of arguments (with NOISETYPE=MGAUSS|MOUTLIERS)");
    } 
  } else sigma_.resize(1, readsigma[0]);

  double read_smin_;
  parse("SIGMA_MIN",read_smin_);
  sigma_min_.resize(sigma_.size(),read_smin_);
  double read_smax_;
  parse("SIGMA_MAX",read_smax_);
  sigma_max_.resize(sigma_.size(),read_smax_);
  double read_dsigma_=-1.;
  parse("DSIGMA",read_dsigma_);
  if(read_dsigma_<0) read_dsigma_ = 0.05*(read_smax_ - read_smin_);
  Dsigma_.resize(sigma_.size(),read_dsigma_);

  // monte carlo stuff
  parse("MC_STEPS",MCsteps_);
  parse("MC_STRIDE",MCstride_);
  parse("MC_CHUNKSIZE", MCchunksize_);
  // adjust for multiple-time steps
  MCstride_ *= getStride();
  // get temperature
  double temp=0.0;
  parse("TEMP",temp);
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();

  // while sigma_mean_ has the same size of sigma
  vector<double> read_sigma_mean_;
  parseVector("SIGMA_MEAN0",read_sigma_mean_);
  if(!do_optsigmamean_ && read_sigma_mean_.size()==0 && !getRestart()) 
    error("If you don't use OPTSIGMAMEAN and you are not RESTARTING then you MUST SET SIGMA_MEAN0");

  if(noise_type_==MGAUSS||noise_type_==MOUTLIERS||noise_type_==GENERIC) {
    if(read_sigma_mean_.size()==narg) {
      sigma_mean2_.resize(narg);
      for(unsigned i=0;i<narg;i++) sigma_mean2_[i]=read_sigma_mean_[i]*read_sigma_mean_[i];
    } else if(read_sigma_mean_.size()==1) {
      sigma_mean2_.resize(narg,read_sigma_mean_[0]*read_sigma_mean_[0]);
    } else if(read_sigma_mean_.size()==0) {
      sigma_mean2_.resize(narg,0.000001);
    } else {
      error("SIGMA_MEAN0 can accept either one single value or as many values as the arguments (with NOISETYPE=MGAUSS|MOUTLIERS)");
    }
  } else {
    if(read_sigma_mean_.size()==1) {
      sigma_mean2_.resize(1, read_sigma_mean_[0]*read_sigma_mean_[0]);
    } else if(read_sigma_mean_.size()==0) {
      sigma_mean2_.resize(narg,0.000001);
    } else {
      error("If you want to use more than one SIGMA_MEAN0 you should use NOISETYPE=MGAUSS|MOUTLIERS");
    }
  } 

  parse("SIGMA_MEAN_CORRECTION", sigma_mean_correction_);

  // sigma mean optimisation
  if(do_optsigmamean_==2) {
    max_force_=3000.;
    parse("MAX_FORCE", max_force_);
    max_force_ *= max_force_;
    sm_mod_=1.0;
    parse("SIGMA_MEAN_MOD0", sm_mod_);
    sm_mod_min_=1.0;
    parse("SIGMA_MEAN_MOD_MIN", sm_mod_min_);
    sm_mod_max_=sqrt(10.);
    parse("SIGMA_MEAN_MOD_MAX", sm_mod_max_);
    Dsm_mod_=0.01;
    parse("DSIGMA_MEAN_MOD", Dsm_mod_);
  }

  checkRead();

  IFile restart_sfile;
  restart_sfile.link(*this);
  if(getRestart()&&restart_sfile.FileExist(status_file_name_)) {
    restart_sfile.open(status_file_name_);
    log.printf("  Restarting from %s\n", status_file_name_.c_str());
    double dummy;
    if(restart_sfile.scanField("time",dummy)){
      for(unsigned i=0;i<sigma_mean2_.size();++i) {
        std::string msg;
        Tools::convert(i,msg);
        double read_sm;
        restart_sfile.scanField("sigma_mean_"+msg,read_sm);
        sigma_mean2_[i]=dummy*dummy;
      }
      for(unsigned i=0;i<sigma_.size();++i) {
        std::string msg;
        Tools::convert(i,msg);
        restart_sfile.scanField("sigma_"+msg,sigma_[i]);
      }
      if(doscale_)  restart_sfile.scanField("scale0_",scale_);
      if(dooffset_) restart_sfile.scanField("offset0_",offset_);
      if(do_optsigmamean_==2) {
        restart_sfile.scanField("sigma_mean_mod0",sm_mod_);
      }
    }
    restart_sfile.scanField();
    restart_sfile.close();

    /* adjust, if needed and wanted, sigma_max  */
    for(unsigned i=0;i<sigma_mean2_.size();++i) {
      double s_v = sqrt(sigma_mean2_[i]*static_cast<double>(nrep_));
      if(do_optsigmamean_>0) {
        if(sigma_max_[i] < s_v) {
          Dsigma_[i] *= s_v/sigma_max_[i]; 
          sigma_max_[i] = s_v;
        }
      }
    }
  }

  switch(noise_type_) {
    case GAUSS:
      log.printf("  with gaussian noise and a single noise parameter for all the data\n");
      break;
    case MGAUSS:
    case GENERIC:
      log.printf("  with gaussian noise and a noise parameter for each data point\n");
      break;
    case OUTLIERS:
      log.printf("  with long tailed gaussian noise and a single noise parameter for all the data\n");
      break;
    case MOUTLIERS:
      log.printf("  with long tailed gaussian noise and a noise parameter for each data point\n");
      break;
  }

  if(doscale_) {
    log.printf("  sampling a common scaling factor with:\n");
    log.printf("    initial scale parameter %f\n",scale_);
    if(scale_prior_==SC_GAUSS) {
      log.printf("    gaussian prior with mean %f and width %f\n",scale_mu_,Dscale_);
    }
    if(scale_prior_==SC_FLAT) {
      log.printf("    flat prior between %f - %f\n",scale_min_,scale_max_);
      log.printf("    maximum MC move of scale parameter %f\n",Dscale_);
    }
  }

  if(dooffset_) {
    log.printf("  sampling a common offset with:\n");
    log.printf("    initial offset parameter %f\n",offset_);
    if(offset_prior_==SC_GAUSS) {
      log.printf("    gaussian prior with mean %f and width %f\n",offset_mu_,Doffset_);
    }
    if(offset_prior_==SC_FLAT) {
      log.printf("    flat prior between %f - %f\n",offset_min_,offset_max_);
      log.printf("    maximum MC move of offset parameter %f\n",Doffset_);
    }
  }

  log.printf("  number of experimental data points %u\n",narg);
  log.printf("  number of replicas %u\n",nrep_);
  log.printf("  initial data uncertainties");
  for(unsigned i=0;i<sigma_.size();++i) log.printf(" %f", sigma_[i]);
  log.printf("\n");
  log.printf("  minimum data uncertainty %f\n",read_smin_);
  log.printf("  maximum data uncertainty %f\n",read_smax_);
  log.printf("  maximum MC move of data uncertainty %f\n",read_dsigma_);
  log.printf("  temperature of the system %f\n",kbt_);
  log.printf("  MC steps %u\n",MCsteps_);
  log.printf("  MC stride %u\n",MCstride_);
  log.printf("  initial standard errors of the mean");
  for(unsigned i=0;i<sigma_mean2_.size();++i) log.printf(" %f", sqrt(sigma_mean2_[i]));
  log.printf("\n");

  if(do_reweight) {
    addComponent("MetaDf");
    componentIsNotPeriodic("MetaDf");
    addComponent("weight");
    componentIsNotPeriodic("weight");
  }

  if(doscale_) { 
    addComponent("scale");  
    componentIsNotPeriodic("scale");
    valueScale=getPntrToComponent("scale");
  }

  if(dooffset_) { 
    addComponent("offset");  
    componentIsNotPeriodic("offset");
    valueOffset=getPntrToComponent("offset");
  }

  if(dooffset_||doscale_) {
    addComponent("acceptScale");
    componentIsNotPeriodic("acceptScale");
    valueAcceptScale=getPntrToComponent("acceptScale");
  }

  addComponent("acceptSigma");
  componentIsNotPeriodic("acceptSigma");
  valueAccept=getPntrToComponent("acceptSigma");

  if(do_optsigmamean_==2) {
    addComponent("maxForceMD");
    componentIsNotPeriodic("maxForceMD");
    valueMaxForceMD=getPntrToComponent("maxForceMD");
    addComponent("smMod");
    componentIsNotPeriodic("smMod");
    valueSMmod=getPntrToComponent("smMod");
  }

  if(noise_type_==MGAUSS||noise_type_==MOUTLIERS||noise_type_==GENERIC) {
    for(unsigned i=0;i<sigma_mean2_.size();++i){
      std::string num; Tools::convert(i,num);
      addComponent("sigmaMean_"+num); componentIsNotPeriodic("sigmaMean_"+num);
      valueSigmaMean.push_back(getPntrToComponent("sigmaMean_"+num));
      getPntrToComponent("sigmaMean_"+num)->set(sqrt(sigma_mean2_[i]));
      addComponent("sigma_"+num); componentIsNotPeriodic("sigma_"+num);
      valueSigma.push_back(getPntrToComponent("sigma_"+num));
      getPntrToComponent("sigma_"+num)->set(sigma_[i]);
      addComponent("ftilde_"+num); componentIsNotPeriodic("ftilde_"+num);
      valueFtilde.push_back(getPntrToComponent("ftilde_"+num));
    }
  } else {
    addComponent("sigmaMean"); componentIsNotPeriodic("sigmaMean");
    valueSigmaMean.push_back(getPntrToComponent("sigmaMean"));
    getPntrToComponent("sigmaMean")->set(sqrt(sigma_mean2_[0]));
    addComponent("sigma"); componentIsNotPeriodic("sigma");
    valueSigma.push_back(getPntrToComponent("sigma"));
    getPntrToComponent("sigma")->set(sigma_[0]);
  }

  // initialize random seed
  unsigned iseed;
  if(master) iseed = time(NULL)+replica_;
  else iseed = 0;     
  comm.Sum(&iseed, 1);
  random[0].setSeed(-iseed);
  // Random chunk
  if(master) iseed = time(NULL)+replica_;
  else iseed = 0;     
  comm.Sum(&iseed, 1);
  random[2].setSeed(-iseed);
  if(doscale_||dooffset_) {
    // in this case we want the same seed everywhere
    iseed = time(NULL);
    if(master) multi_sim_comm.Bcast(iseed,0);
    comm.Bcast(iseed,0);
    random[1].setSeed(-iseed);
  }

  // outfile stuff
  if(write_stride_>0) {
    sfile_.link(*this);
    sfile_.open(status_file_name_);
  }

  log<<"  Bibliography "<<plumed.cite("Bonomi, Camilloni, Cavalli, Vendruscolo, Sci. Adv. 2, e150117 (2016)");
  if(do_reweight) log<<plumed.cite("Bonomi, Camilloni, Vendruscolo, Sci. Rep. 6, 31232 (2016)");
  log<<"\n";
}

Metainference::~Metainference()
{
  if(sfile_.isOpen()) sfile_.close();
}

double Metainference::getEnergySP(const vector<double> &mean, const vector<double> &sigma, 
                                  const double scale, const double offset, const double modifier)
{
  const double scale2 = scale*scale;
  const double mod2   = modifier*modifier;
  const double sm2    = mod2*sigma_mean2_[0];
  const double ss2    = sigma[0]*sigma[0] + scale2*sm2;
  const double sss    = sigma[0]*sigma[0] + sm2;

  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  {
    #pragma omp for reduction( + : ene)
    for(unsigned i=0;i<narg;++i){
      const double dev = scale*mean[i]-parameters[i]+offset; 
      const double a2 = 0.5*dev*dev + ss2;
      ene += std::log(2.0*a2/(1.0-exp(-a2/sm2)));
    }
  }
  // add one single Jeffrey's prior and one normalisation per data point
  ene += std::log(sss) + static_cast<double>(narg)*0.5*std::log(0.5*M_PI*M_PI/ss2);
  if(doscale_)  ene += std::log(sss);
  if(dooffset_) ene += std::log(sss);
  return kbt_ * ene;
}

double Metainference::getEnergySPE(const vector<double> &mean, const vector<double> &sigma, 
                                   const double scale, const double offset, const double modifier)
{
  const double scale2 = scale*scale;
  const double mod2   = modifier*modifier;
  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  {
    #pragma omp for reduction( + : ene)
    for(unsigned i=0;i<narg;++i){
      const double sm2 = mod2*sigma_mean2_[i];
      const double ss2 = sigma[i]*sigma[i] + scale2*sm2;
      const double sss = sigma[i]*sigma[i] + sm2;
      const double dev = scale*mean[i]-parameters[i]+offset; 
      const double a2  = 0.5*dev*dev + ss2;
      ene += std::log(sss) + 0.5*std::log(0.5*M_PI*M_PI/ss2) + std::log(2.0*a2/(1.0-exp(-a2/sm2)));
      if(doscale_)  ene += 0.5*std::log(sss);
      if(dooffset_) ene += 0.5*std::log(sss);
    }
  }
  return kbt_ * ene;
}

double Metainference::getEnergyMGGEN(const vector<double> &mean, const vector<double> &ftilde, const vector<double> &sigma, 
                                     const double scale, const double offset, const double modifier)
{
  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  { 
    #pragma omp for reduction( + : ene)
    for(unsigned i=0;i<narg;++i){
      const double inv_sb2  = 1./(sigma[i]*sigma[i]);
      const double inv_sm2  = 1./(sigma_mean2_[i]);
      double devb = scale*ftilde[i]-parameters[i]+offset;
      double devm = mean[i] - ftilde[i];
      // deviation + normalisation + jeffrey
      const double normb         = -0.5*std::log(0.5/M_PI*inv_sb2);
      const double normm         = -0.5*std::log(0.5/M_PI*inv_sm2);
      const double jeffreys      = -0.5*std::log(2.*inv_sb2);
      ene += 0.5*devb*devb*inv_sb2 + 0.5*devm*devm*inv_sm2 + normb + normm + jeffreys;
      if(doscale_)  ene += jeffreys;
      if(dooffset_) ene += jeffreys;
    }
  }
  return kbt_ * ene;
}

double Metainference::getEnergyGJ(const vector<double> &mean, const vector<double> &sigma, 
                                  const double scale, const double offset, const double modifier)
{
  const double scale2  = scale*scale;
  const double mod2    = modifier*modifier;
  const double inv_s2  = 1./(sigma[0]*sigma[0] + scale2*mod2*sigma_mean2_[0]);
  const double inv_sss = 1./(sigma[0]*sigma[0] + mod2*sigma_mean2_[0]);

  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  {
    #pragma omp for reduction( + : ene)
    for(unsigned i=0;i<narg;++i){
      double dev = scale*mean[i]-parameters[i]+offset;
      ene += 0.5*dev*dev*inv_s2;
    }
  }
  const double normalisation = -0.5*std::log(0.5/M_PI*inv_s2);
  const double jeffreys = -0.5*std::log(2.*inv_sss);
  // add Jeffrey's prior in case one sigma for all data points + one normalisation per datapoint
  ene += jeffreys + static_cast<double>(narg)*normalisation;
  if(doscale_)  ene += jeffreys;
  if(dooffset_) ene += jeffreys;

  return kbt_ * ene;
}

double Metainference::getEnergyGJE(const vector<double> &mean, const vector<double> &sigma, 
                                   const double scale, const double offset, const double modifier)
{
  const double scale2 = scale*scale;
  const double mod2   = modifier*modifier;

  double ene = 0.0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene)
  { 
    #pragma omp for reduction( + : ene)
    for(unsigned i=0;i<narg;++i){
      const double inv_s2  = 1./(sigma[i]*sigma[i] + scale2*mod2*sigma_mean2_[i]);
      const double inv_sss = 1./(sigma[i]*sigma[i] + mod2*sigma_mean2_[i]);
      double dev = scale*mean[i]-parameters[i]+offset;
      // deviation + normalisation + jeffrey
      const double normalisation = -0.5*std::log(0.5/M_PI*inv_s2);
      const double jeffreys      = -0.5*std::log(2.*inv_sss);
      ene += 0.5*dev*dev*inv_s2 + normalisation + jeffreys;
      if(doscale_)  ene += jeffreys;
      if(dooffset_) ene += jeffreys;
    }
  }
  return kbt_ * ene;
}

void Metainference::doMonteCarlo(const vector<double> &mean_, const double modifier)
{
  // calculate old energy with the updated coordinates
  double old_energy=0.;

  switch(noise_type_) {
    case GAUSS:
      old_energy = getEnergyGJ(mean_,sigma_,scale_,offset_,modifier);
      break;
    case MGAUSS:
      old_energy = getEnergyGJE(mean_,sigma_,scale_,offset_,modifier);
      break;
    case OUTLIERS:
      old_energy = getEnergySP(mean_,sigma_,scale_,offset_,modifier);
      break;
    case MOUTLIERS:
      old_energy = getEnergySPE(mean_,sigma_,scale_,offset_,modifier);
      break;
    case GENERIC:
      old_energy = getEnergyMGGEN(mean_,ftilde_,sigma_,scale_,offset_,modifier);
      break;
  }

  // Create vector of random sigma indices
  vector<unsigned> indices;
  if (MCchunksize_ > 0) {
      for (unsigned j=0; j<sigma_.size(); j++) {
          indices.push_back(j);
      }
      random[2].Shuffle(indices);
  }
  bool breaknow = false;

  // cycle on MC steps 
  for(unsigned i=0;i<MCsteps_;++i){
 
    MCtrial_++;

    // propose move for ftilde 
    vector<double> new_ftilde(sigma_.size());
    new_ftilde = ftilde_;

    if(noise_type_==GENERIC) {
      // change all sigmas
      for(unsigned j=0;j<sigma_.size();j++) {
        const double r3 = random[0].Gaussian();
        //const double ds3 = 0.1*sqrt(sigma_mean_[j]*sigma_mean_[j]+sigma_[j]*sigma_[j])*r3;
        const double ds3 = 0.2*sqrt(sigma_mean2_[j])*r3;
        new_ftilde[j] = ftilde_[j] + ds3;
      }
      // calculate new energy
      double new_energy = getEnergyMGGEN(mean_,new_ftilde,sigma_,scale_,offset_,modifier);

      // accept or reject
      const double delta = ( new_energy - old_energy ) / kbt_;
      // if delta is negative always accept move
      if( delta <= 0.0 ){
        old_energy = new_energy;
        ftilde_ = new_ftilde;
        MCacceptFT_++;
      // otherwise extract random number
      } else {
        const double s = random[0].RandU01();
        if( s < exp(-delta) ){
          old_energy = new_energy;
          ftilde_ = new_ftilde;
          MCacceptFT_++;
        }
      }
    }

    // propose move for scale and/or offset
    double new_scale = scale_;
    double new_offset = offset_;

    if(doscale_||dooffset_) {
      if(doscale_) {
        if(scale_prior_==SC_FLAT) {
          const double r1 = random[1].Gaussian();
          const double ds1 = Dscale_*r1;
          new_scale += ds1;
          // check boundaries
          if(new_scale > scale_max_){new_scale = 2.0 * scale_max_ - new_scale;}
          if(new_scale < scale_min_){new_scale = 2.0 * scale_min_ - new_scale;}
        } else {
          const double r1 = random[1].Gaussian();
          const double ds1 = 0.5*(scale_mu_-new_scale)+Dscale_*exp(1)/M_PI*r1;
          new_scale += ds1;
        }
      }

      if(dooffset_) {
        if(offset_prior_==SC_FLAT) {
          const double r1 = random[1].Gaussian();
          const double ds1 = Doffset_*r1;
          new_offset += ds1;
          // check boundaries
          if(new_offset > offset_max_){new_offset = 2.0 * offset_max_ - new_offset;}
          if(new_offset < offset_min_){new_offset = 2.0 * offset_min_ - new_offset;}
        } else {
          const double r1 = random[1].Gaussian();
          const double ds1 = 0.5*(offset_mu_-new_offset)+Doffset_*exp(1)/M_PI*r1;
          new_offset += ds1;
        }
      }

      // calculate new energy
      double new_energy = 0.;

      switch(noise_type_) {
        case GAUSS:
          new_energy = getEnergyGJ(mean_,sigma_,new_scale,new_offset,modifier);
          break;
        case MGAUSS:
          new_energy = getEnergyGJE(mean_,sigma_,new_scale,new_offset,modifier);
          break;
        case OUTLIERS:
          new_energy = getEnergySP(mean_,sigma_,new_scale,new_offset,modifier);
          break;
        case MOUTLIERS:
          new_energy = getEnergySPE(mean_,sigma_,new_scale,new_offset,modifier);
          break;
        case GENERIC:
         new_energy = getEnergyMGGEN(mean_,ftilde_,sigma_,new_scale,new_offset,modifier);
         break;
      }
      // for the scale we need to consider the total energy
      vector<double> totenergies(2);
      if(master) {
        totenergies[0] = old_energy;
        totenergies[1] = new_energy;
        if(nrep_>1) multi_sim_comm.Sum(totenergies);
      } else {
        totenergies[0] = 0;
        totenergies[1] = 0;
      }
      comm.Sum(totenergies);

      // accept or reject
      const double delta = ( totenergies[1] - totenergies[0] ) / kbt_;
      // if delta is negative always accept move
      if( delta <= 0.0 ){
        old_energy = new_energy;
        scale_ = new_scale;
        offset_ = new_offset;
        MCacceptScale_++;
      // otherwise extract random number
      } else {
        double s = random[1].RandU01();
        if( s < exp(-delta) ){
          old_energy = new_energy;
          scale_ = new_scale;
          offset_ = new_offset;
          MCacceptScale_++;
        }
      }
    }

    // propose move for sigma
    vector<double> new_sigma(sigma_.size());
    new_sigma = sigma_;

    // change MCchunksize_ sigmas
    if (MCchunksize_ > 0) {
        if ((MCchunksize_ * i) >= sigma_.size()) {
            // This means we are not moving any sigma, so we should break immediately
            breaknow = true;
        }

        // change random sigmas
        for(unsigned j=0;j<MCchunksize_;j++) {
            const unsigned shuffle_index = j + MCchunksize_ * i;
            if (shuffle_index >= sigma_.size()) {
                // Going any further will segfault but we should still evaluate the sigmas we changed
                break;
            }
            const unsigned index = indices[shuffle_index];
            const double r2 = random[0].Gaussian();
            const double ds2 = Dsigma_[index]*r2;
            new_sigma[index] = sigma_[index] + ds2;
            // check boundaries
            if(new_sigma[index] > sigma_max_[index]){new_sigma[index] = 2.0 * sigma_max_[index] - new_sigma[index];}
            if(new_sigma[index] < sigma_min_[index]){new_sigma[index] = 2.0 * sigma_min_[index] - new_sigma[index];}
        }
    } else {
        // change all sigmas
        for(unsigned j=0;j<sigma_.size();j++) {
            const double r2 = random[0].Gaussian();
            const double ds2 = Dsigma_[j]*r2;
            new_sigma[j] = sigma_[j] + ds2;
            // check boundaries
            if(new_sigma[j] > sigma_max_[j]){new_sigma[j] = 2.0 * sigma_max_[j] - new_sigma[j];}
            if(new_sigma[j] < sigma_min_[j]){new_sigma[j] = 2.0 * sigma_min_[j] - new_sigma[j];}
        }
    }

    if (breaknow) {
        // We didnt move any sigmas, so no sense in evaluating anything
        break;
    }

    // calculate new energy
    double new_energy;
    switch(noise_type_) {
      case GAUSS:
        new_energy = getEnergyGJ(mean_,new_sigma,scale_,offset_,modifier);
        break;
      case MGAUSS:
        new_energy = getEnergyGJE(mean_,new_sigma,scale_,offset_,modifier);
        break;
      case OUTLIERS:
        new_energy = getEnergySP(mean_,new_sigma,scale_,offset_,modifier);
        break;
      case MOUTLIERS:
        new_energy = getEnergySPE(mean_,new_sigma,scale_,offset_,modifier);
        break;
      case GENERIC:
        new_energy = getEnergyMGGEN(mean_,ftilde_,new_sigma,scale_,offset_,modifier);
        break;
    }

    // accept or reject
    const double delta = ( new_energy - old_energy ) / kbt_;
    // if delta is negative always accept move
    if( delta <= 0.0 ){
      old_energy = new_energy;
      sigma_ = new_sigma;
      MCaccept_++;
    // otherwise extract random number
    } else {
      const double s = random[0].RandU01();
      if( s < exp(-delta) ){
        old_energy = new_energy;
        sigma_ = new_sigma;
        MCaccept_++;
      }
    }
 
  }
  /* save the result of the sampling */
  double accept = static_cast<double>(MCaccept_) / static_cast<double>(MCtrial_);
  valueAccept->set(accept);
  if(doscale_) valueScale->set(scale_);
  if(dooffset_) valueOffset->set(offset_);
  if(doscale_||dooffset_) {
    accept = static_cast<double>(MCacceptScale_) / static_cast<double>(MCtrial_);
    valueAcceptScale->set(accept);
  }
  for(unsigned i=0; i<sigma_.size(); i++) valueSigma[i]->set(sigma_[i]);
  for(unsigned i=0; i<sigma_.size(); i++) valueFtilde[i]->set(ftilde_[i]);
}

/* 
   In the following energy-force functions we don't add the normalisation and the jeffreys priors
   because they are not needed for the forces, the correct MetaInference energy is the one calculated
   in the Monte-Carlo 
*/

double Metainference::getEnergyForceSP(const vector<double> &mean, const vector<double> &dmean_x,
                                       const vector<double> &dmean_b, const double modifier)
{
  const double mod2   = modifier*modifier; /* this is now modifiers */
  const double scale2 = scale_*scale_;
  const double sm2    = mod2*sigma_mean2_[0]; 
  const double ss2    = sigma_[0]*sigma_[0] + scale2*sm2;
  vector<double> f(narg+1,0);
  
  if(master){
    double omp_ene=0.;
    #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(omp_ene)
    { 
      #pragma omp for reduction( + : omp_ene)
      for(unsigned i=0;i<narg;++i){
        const double dev = scale_*mean[i]-parameters[i]+offset_; 
        const double a2 = 0.5*dev*dev + ss2;
        const double t = exp(-a2/sm2);
        const double dt = 1./t;
        const double it = 1./(1.-t);
        const double dit = 1./(1.-dt);
        omp_ene += std::log(2.*a2*it);
        f[i] = -scale_*dev*(dit/sm2 + 1./a2);
      }
    }
    f[narg] = omp_ene;
    // collect contribution to forces and energy from other replicas
    if(nrep_>1) multi_sim_comm.Sum(&f[0],narg+1);
  }
  // intra-replica summation
  comm.Sum(&f[0],narg+1);

  const double ene = f[narg];
  double w_tmp = 0.;
  for(unsigned i=0; i<narg; ++i) {
    setOutputForce(i, kbt_*f[i]*dmean_x[i]);
    w_tmp += kbt_*f[i]*dmean_b[i];
  }

  if(do_reweight) {
    setOutputForce(narg, w_tmp);
    getPntrToComponent("MetaDf")->set(-w_tmp);
  }

  return kbt_*ene;
}

double Metainference::getEnergyForceSPE(const vector<double> &mean, const vector<double> &dmean_x,
                                        const vector<double> &dmean_b, const double modifier)
{
  const double mod2   = modifier*modifier;
  const double scale2 = scale_*scale_;
  vector<double> f(narg+1,0);

  if(master){
    double omp_ene = 0;
    #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(omp_ene)
    { 
      #pragma omp for reduction( + : omp_ene)
      for(unsigned i=0;i<narg;++i){
        const double sm2 = mod2*sigma_mean2_[i]; 
        const double ss2 = sigma_[i]*sigma_[i] + scale2*sm2;
        const double dev = scale_*mean[i]-parameters[i]+offset_; 
        const double a2  = 0.5*dev*dev + ss2;
        const double t   = exp(-a2/sm2);
        const double dt  = 1./t;
        const double it  = 1./(1.-t);
        const double dit = 1./(1.-dt);
        omp_ene += std::log(2.*a2*it);
        f[i] = -scale_*dev*(dit/sm2 + 1./a2);
      }
    }
    f[narg] = omp_ene;
    // collect contribution to forces and energy from other replicas
    if(nrep_>1) multi_sim_comm.Sum(&f[0],narg+1);
  }
  comm.Sum(&f[0],narg+1);

  const double ene = f[narg];
  double w_tmp = 0.;
  for(unsigned i=0; i<narg; ++i) {
    setOutputForce(i, kbt_ * dmean_x[i] * f[i]);
    w_tmp += kbt_ * dmean_b[i] *f[i];
  }

  if(do_reweight) {
    setOutputForce(narg, w_tmp);
    getPntrToComponent("MetaDf")->set(-w_tmp);
  }

  return kbt_*ene;
}

double Metainference::getEnergyForceGJ(const vector<double> &mean, const vector<double> &dmean_x,
                                       const vector<double> &dmean_b, const double modifier)
{
  const double mod2   = modifier*modifier;
  const double scale2 = scale_*scale_;
  double inv_s2=0.;

  if(master) {
    inv_s2 = 1./(sigma_[0]*sigma_[0] + mod2*scale2*sigma_mean2_[0]);
    if(nrep_>1) multi_sim_comm.Sum(inv_s2);
  } 
  comm.Sum(inv_s2);  

  double ene   = 0.;
  double w_tmp = 0.;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene,w_tmp)
  { 
    #pragma omp for reduction( + : ene,w_tmp)
    for(unsigned i=0;i<narg;++i){
      const double dev = scale_*mean[i]-parameters[i]+offset_;
      const double mult = dev*scale_*inv_s2;
      ene += 0.5*dev*dev*inv_s2;
      setOutputForce(i, -kbt_*dmean_x[i]*mult);
      w_tmp += kbt_*dmean_b[i]*mult;
    }
  }

  if(do_reweight) {
    setOutputForce(narg, -w_tmp);
    getPntrToComponent("MetaDf")->set(w_tmp);
  }

  return kbt_*ene;
}

double Metainference::getEnergyForceGJE(const vector<double> &mean, const vector<double> &dmean_x, const vector<double> &dmean_b, 
                                        const vector<double> &dsigma_mean2_x, const vector<double> &dsigma_mean2_b, const double modifier)
{
  const double   mod2   = modifier*modifier;
  const double   scale2 = scale_*scale_;
  vector<double> inv_s2(sigma_.size(),0.);
  vector<double> inv_ss(sigma_.size(),0.);
  vector<double> inv2_s2(sigma_.size(),0.);
  vector<double> normg(sigma_.size(),0.);
  vector<double> jeff(sigma_.size(),0.);

  if(master) {
    for(unsigned i=0;i<sigma_.size(); ++i) {
      inv_s2[i]  = 1./(sigma_[i]*sigma_[i] + scale2*mod2*sigma_mean2_[i]);
      inv_ss[i]  = 1./(sigma_[i]*sigma_[i] + mod2*sigma_mean2_[i]);
      inv2_s2[i] = inv_s2[i]*inv_s2[i];
      normg[i] = -0.5*std::log(0.5/M_PI*inv_s2[i]);
      jeff[i]  = -0.5*std::log(2.*inv_ss[i]);
    }
    if(nrep_>1) {
      multi_sim_comm.Sum(&inv_s2[0],sigma_.size());
      multi_sim_comm.Sum(&inv_ss[0],sigma_.size());
      multi_sim_comm.Sum(&inv2_s2[0],sigma_.size());
      multi_sim_comm.Sum(&normg[0],sigma_.size());
      multi_sim_comm.Sum(&jeff[0],sigma_.size());
    }
  }
  comm.Sum(&inv_s2[0],sigma_.size());  
  comm.Sum(&inv_ss[0],sigma_.size());  
  comm.Sum(&inv2_s2[0],sigma_.size());  
  comm.Sum(&normg[0],sigma_.size());
  comm.Sum(&jeff[0],sigma_.size());
 
  double dene_b = 0.;
  double ene    = 0.;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene,dene_b)
  { 
    #pragma omp for reduction( + : ene,dene_b)
    for(unsigned i=0;i<narg;++i){
      const double dev = scale_*mean[i]-parameters[i]+offset_;
      double dene_x  = (kbt_*scale_*dev*inv_s2[i]*dmean_x[i] - 0.5*kbt_*dev*dev*inv2_s2[i]*scale2*mod2*dsigma_mean2_x[i]);
             dene_b += (kbt_*scale_*dev*inv_s2[i]*dmean_b[i] - 0.5*kbt_*dev*dev*inv2_s2[i]*scale2*mod2*dsigma_mean2_b[i]);
             dene_x += 0.5*kbt_*scale2*mod2*inv_s2[i]*dsigma_mean2_x[i];
             dene_b += 0.5*kbt_*scale2*mod2*inv_s2[i]*dsigma_mean2_b[i];
             dene_x += 0.5*kbt_*mod2*inv_ss[i]*dsigma_mean2_x[i];
             dene_b += 0.5*kbt_*mod2*inv_ss[i]*dsigma_mean2_b[i];
      ene += 0.5*dev*dev*inv_s2[i] + normg[i] + jeff[i];
      setOutputForce(i, -dene_x);
    }
  }

  if(do_reweight) {
    setOutputForce(narg, -dene_b);
    getPntrToComponent("MetaDf")->set(dene_b);
  }

  return kbt_*ene;
}

double Metainference::getEnergyForceMGGEN(const vector<double> &mean, const vector<double> &dmean_x, const vector<double> &dmean_b, 
                                          const vector<double> &dsigma_mean2_x, const vector<double> &dsigma_mean2_b, const double modifier)
{
  vector<double> inv_s2(sigma_.size(),0.);
  vector<double> dev(sigma_.size(),0.);
  vector<double> dev2(sigma_.size(),0.);
  vector<double> normg(sigma_.size(),0.);

  for(unsigned i=0;i<sigma_.size(); ++i) {
    inv_s2[i]   = 1./(sigma_mean2_[i]);
    normg[i]    = -0.5*nrep_*std::log(0.5/M_PI*inv_s2[i]);
    if(master) {
      dev[i]  = (mean[i]-ftilde_[i]);
      dev2[i] = (mean[i]-ftilde_[i])*(mean[i]-ftilde_[i]);
    }
  }
  if(master&&nrep_>1) {
    multi_sim_comm.Sum(&dev[0],dev.size());
    multi_sim_comm.Sum(&dev2[0],dev.size());
  }
  comm.Sum(&dev[0],dev.size());
  comm.Sum(&dev2[0],dev.size());
 
  double dene_b = 0.;
  double ene    = 0.;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(ene,dene_b)
  { 
    #pragma omp for reduction( + : ene,dene_b)
    for(unsigned i=0;i<narg;++i){
      double dene_x  = kbt_*inv_s2[i]*dmean_x[i]*dev[i] - 0.5*kbt_*inv_s2[i]*inv_s2[i]*dev2[i]*dsigma_mean2_x[i];
             dene_b += kbt_*inv_s2[i]*dmean_b[i]*dev[i] - 0.5*kbt_*inv_s2[i]*inv_s2[i]*dev2[i]*dsigma_mean2_b[i];
             dene_x += 0.5*kbt_*nrep_*inv_s2[i]*dsigma_mean2_x[i];
             dene_b += 0.5*kbt_*nrep_*inv_s2[i]*dsigma_mean2_b[i];
      ene += 0.5*dev2[i]*inv_s2[i] + normg[i];
      setOutputForce(i, -dene_x);
    }
  }

  if(do_reweight) {
    setOutputForce(narg, -dene_b);
    getPntrToComponent("MetaDf")->set(dene_b);
  }

  return kbt_*ene;
}

void Metainference::calculate()
{
  const double dnrep = static_cast<double>(nrep_);
  double norm        = 0.0;
  double vnorm       = 0.0;
  double fact        = 0.0;
  vector<double> sigma_mean2(narg,0);
  vector<double> dsigma_mean2_x(narg,0);
  vector<double> dsigma_mean2_b(narg,0);

  if(do_reweight){
    // calculate the weights either from BIAS 
    vector<double> bias(nrep_,0);
    if(master){
      bias[replica_] = getArgument(narg); 
      if(nrep_>1) multi_sim_comm.Sum(&bias[0], nrep_);  
    }
    comm.Sum(&bias[0], nrep_);
    const double maxbias = *(std::max_element(bias.begin(), bias.end()));
    for(unsigned i=0; i<nrep_; ++i){
      bias[i] = exp((bias[i]-maxbias)/kbt_); 
      norm   += bias[i];
    }
    for(unsigned i=0;i<nrep_;++i) vnorm += (bias[i]/norm)*(bias[i]/norm);
    fact = bias[replica_]/norm;
    getPntrToComponent("weight")->set(fact);
  } else {
    // or arithmetic ones
    norm = dnrep; 
    fact = 1.0/norm; 
  }

  // calculate the mean 
  vector<double> mean(narg,0);
  // this is the derivative of the mean with respect to the argument
  vector<double> dmean_x(narg,fact);
  // this is the derivative of the mean with respect to the bias
  vector<double> dmean_b(narg,0);
  if(master) {
    for(unsigned i=0;i<narg;++i) mean[i] = fact*getArgument(i); 
    if(nrep_>1) multi_sim_comm.Sum(&mean[0], narg);
  }
  comm.Sum(&mean[0], narg);
  // set the derivative of the mean with respect to the bias
  for(unsigned i=0;i<narg;++i) dmean_b[i] = fact/kbt_*(getArgument(i)-mean[i]);

  if(firsttime) {ftilde_ = mean; firsttime = false;}

  for(unsigned i=0;i<narg;++i) dmean_b[i] = fact/kbt_*(getArgument(i)-mean[i]);

  if(do_optsigmamean_>0) {
    vector<double> v_tmp1(narg,0);
    vector<double> v_tmp2(narg,0);
    /* this is the current estimate of sigma mean for each argument
       there is one of this per argument in any case  because it is
       the maximum among these to be used in case of GAUSS/OUTLIER */
    if(do_reweight) {
      if(master) {
        for(unsigned i=0;i<narg;++i) { 
          v_tmp1[i] = fact*(getArgument(i)-mean[i])*(getArgument(i)-mean[i]);
          v_tmp2[i] = fact*(getArgument(i)-mean[i]);
        }
        if(nrep_>1) {
          multi_sim_comm.Sum(&v_tmp1[0], narg);
          multi_sim_comm.Sum(&v_tmp2[0], narg);
        }
      }
      comm.Sum(&v_tmp1[0], narg);
      comm.Sum(&v_tmp2[0], narg);
      for(unsigned i=0;i<narg;++i) {
        sigma_mean2[i] = 2.*vnorm/(1.-vnorm)*v_tmp1[i];
        dsigma_mean2_x[i] = 2.*2.*vnorm/(1.-vnorm)*(-v_tmp2[i]*dmean_x[i]+fact*(getArgument(i)-mean[i]));
        double part1 = +2.*vnorm/((1.-vnorm)*(1.-vnorm)*kbt_)*(-fact*fact+fact*vnorm)*v_tmp1[i];
        double part2 = +2./((1.-vnorm)*kbt_)*fact*fact*v_tmp1[i];
        double part3 = -2.*vnorm*fact/((1.-vnorm)*kbt_)*v_tmp1[i];
        double part4 = +2.*vnorm/((1.-vnorm))*(-fact*v_tmp1[i]/kbt_+fact/kbt_*(getArgument(i)-mean[i])*(getArgument(i)-mean[i])-2.*dmean_b[i]*v_tmp2[i]);
        dsigma_mean2_b[i] = 2.*(part1+part2+part3+part4);
      }
    } else {
      /* standard estimate of sigma_mean currently missing */
    }

    if(noise_type_==MGAUSS||noise_type_==MOUTLIERS) {
      for(unsigned i=0;i<narg;++i) {
        /* the standard error of the mean */
        valueSigmaMean[i]->set(sqrt(sigma_mean2_[i]));
        /* this is the variance */
        const double s_v = sqrt(sigma_mean2_[i]*dnrep);
        /* if sigma_max is less than the variance we increase it and increase Dsigma accordingly */
        if(sigma_max_[i] < s_v) {
          Dsigma_[i] *= s_v/sigma_max_[i]; 
          sigma_max_[i] = s_v;
        }
        sigma_min_[i] = sqrt(sigma_mean2_[i]);
        if(sigma_[i] < sigma_min_[i]) sigma_[i] = sigma_min_[i];
      }
    }
  }

  // rescale sigma_mean by a user supplied constant 
  double sigma_mean_modifier = sigma_mean_correction_;
  /* fix sigma_mean_ for the effect of large forces */
  if(do_optsigmamean_==2) sigma_mean_modifier *= sm_mod_;

  /* MONTE CARLO */
  const long int step = getStep();
  if(step%MCstride_==0&&!getExchangeStep()) doMonteCarlo(mean, sigma_mean_modifier);

  // calculate bias and forces
  double ene = 0; 
  switch(noise_type_) {
    case GAUSS:
      ene = getEnergyForceGJ(mean, dmean_x, dmean_b, sigma_mean_modifier);
      break;
    case MGAUSS:
      ene = getEnergyForceGJE(mean, dmean_x, dmean_b, dsigma_mean2_x, dsigma_mean2_b, sigma_mean_modifier);
      break;
    case OUTLIERS:
      ene = getEnergyForceSP(mean, dmean_x, dmean_b, sigma_mean_modifier);
      break;
    case MOUTLIERS:
      ene = getEnergyForceSPE(mean, dmean_x, dmean_b, sigma_mean_modifier);
      break;
    case GENERIC:
      ene = getEnergyForceMGGEN(mean, dmean_x, dmean_b, dsigma_mean2_x, dsigma_mean2_b, sigma_mean_modifier);
      break;
  }

  // set value of the bias
  setBias(ene);
}

void Metainference::writeStatus()
{
  sfile_.rewind();
  sfile_.printField("time",getTimeStep()*getStep());
  for(unsigned i=0;i<sigma_mean2_.size();++i) {
    std::string msg;
    Tools::convert(i,msg);
    sfile_.printField("sigma_mean_"+msg,sqrt(sigma_mean2_[i]));
  }
  for(unsigned i=0;i<sigma_.size();++i) {
    std::string msg;
    Tools::convert(i,msg);
    sfile_.printField("sigma_"+msg,sigma_[i]);
  }
  if(doscale_) {
    sfile_.printField("scale0_",scale_);
  }
  if(dooffset_) {
    sfile_.printField("offset0_",offset_);
  }
  if(do_optsigmamean_==2) {
    sfile_.printField("sigma_mean_mod0",sm_mod_);
  }
  sfile_.printField();
  sfile_.flush();
}

void Metainference::update() {
  if(do_optsigmamean_==2) {
    const double EPS = 0.1;
    // Get max force of whole system
    vector<Vector> md_forces;
    vector<Vector> plumed_forces;
    vector<double> masses;

    atoms.getLocalMDForces(md_forces);
    atoms.getLocalForces(plumed_forces);
    atoms.getLocalMasses(masses);

    vector<double> allforces_md;
    allforces_md.reserve(md_forces.size());

    for(unsigned i = 0; i < plumed_forces.size(); ++i) {
      const double pf2 = plumed_forces[i].modulo2();
      // we are only interested in forces plumed has an effect on
      if(pf2 > EPS && masses[i] > EPS ) {
        const double invm2 = 1./(masses[i]*masses[i]);
        allforces_md.push_back(md_forces[i].modulo2()*invm2);
      }
    }

    vector<double> fmax_tmp_md(comm.Get_size(),0);
    vector<double> nfmax_md(nrep_, 0);

    /* each local thread should look for the maximum force and number of violations */
    if(allforces_md.size()>0) {
      fmax_tmp_md[comm.Get_rank()] = *max_element(allforces_md.begin(), allforces_md.end());
      for(unsigned i = 0; i < allforces_md.size(); ++i) {
        if(allforces_md[i] > max_force_) {
          nfmax_md[replica_] += allforces_md[i]/max_force_;
        }
      }
    }
    // the largest forces are shared among the local threads but not over the replicas 
    comm.Sum(fmax_tmp_md);
    // these are the largest forces for a specific replica
    const double fmax_md = *max_element(fmax_tmp_md.begin(), fmax_tmp_md.end());

    // the number of violations is summed up over the local thread and over the replicas 
    comm.Sum(nfmax_md);
    if(master && nrep_ > 1) multi_sim_comm.Sum(nfmax_md);
    comm.Bcast(&nfmax_md[0], nrep_, 0);
    
    const double nnfmax_md = (*max_element(nfmax_md.begin(), nfmax_md.end()))*nrep_;

    if( nnfmax_md == 0) {
      sm_mod_ -= Dsm_mod_ * 0.01 * std::log(sm_mod_/sm_mod_min_);
      if(sm_mod_<sm_mod_min_) sm_mod_=sm_mod_min_;
    } else {
      const double sm_mod_new = sm_mod_ + Dsm_mod_ * std::log(nnfmax_md+1.);
      if(sm_mod_new > sm_mod_max_) {
        sm_mod_ = sm_mod_max_;
      } else {
        sm_mod_ = sm_mod_new;
      }
    }

    valueSMmod->set(sm_mod_);
    valueMaxForceMD->set(sqrt(fmax_md));
  }

  // write status file
  const long int step = getStep();
  if(write_stride_>0&& (step%write_stride_==0 || getCPT()) ) writeStatus();
}

}
}

