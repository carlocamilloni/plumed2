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

*/
#include "bias/Bias.h"
#include "bias/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "core/Value.h"
#include "tools/File.h"
#include "tools/Random.h"
#include <cmath>
#include <ctime>

using namespace std;

namespace PLMD{
namespace isdb{

//+PLUMEDOC BIAS RESCALE
/*

*/
//+ENDPLUMEDOC

class Rescale : public bias::Bias
{
  // beta parameter
  vector<double> beta_;
  unsigned       ibeta_;
  double         w0_;
  double         biasf_;
  vector<double> bias_;
  vector<double> expo_;
  vector<unsigned> shared_;
  // print bias
  unsigned int   Biasstride_;
  string         Biasfilename_;
  bool           first_bias_;
  OFile          Biasfile_;
  // temperature in kbt
  double kbt_;
  // Monte Carlo stuff
  unsigned MCsteps_;
  unsigned MCstride_;
  long int MCfirst_;
  long unsigned MCaccbeta_;
  // replica stuff
  unsigned nrep_;
  unsigned replica_;

  // Monte Carlo
  void doMonteCarlo(long int step, vector<double> args);
  unsigned  proposeMove(unsigned x, unsigned xmin, unsigned xmax);
  bool doAccept(double oldE, double newE, double oldB, double newB);
  // read and print bias
  void read_bias();
  void print_bias(long int step);
  
public:
  explicit Rescale(const ActionOptions&);
  ~Rescale();
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Rescale,"RESCALE")

void Rescale::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","TEMP","temperature");
  keys.add("compulsory","IBETA0","initial value of the ibeta parameter");
  keys.add("compulsory","BETA_MAX","maximum values of the beta parameter");
  keys.add("compulsory","NBIN","number of bins for beta grid");
  keys.add("compulsory","W0", "initial bias height");
  keys.add("compulsory","BIASFACTOR", "bias factor");
  keys.add("compulsory","BSTRIDE", "stride for writing bias");
  keys.add("compulsory","BFILE", "file name for bias");
  keys.add("compulsory","SHARED", "mark arguments that need to be summed across replicas");
  keys.add("optional","MC_STEPS","number of MC steps");
  keys.add("optional","MC_STRIDE","MC stride");
  keys.addOutputComponent("ibeta",   "default","ibeta parameter");
  keys.addOutputComponent("accbeta", "default","MC acceptance beta");
  keys.addOutputComponent("wtbias",  "default","well-tempered bias");
}

Rescale::Rescale(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao), 
first_bias_(true),
MCsteps_(1), MCstride_(1), MCfirst_(-1), MCaccbeta_(0)
{
  // set up replica stuff 
  if(comm.Get_rank()==0) {
    nrep_    = multi_sim_comm.Get_size();
    replica_ = multi_sim_comm.Get_rank();
  } else {
    nrep_    = 0;
    replica_ = 0;
  }
  comm.Sum(&nrep_,1);
  comm.Sum(&replica_,1);

  // beta parameter
  parse("IBETA0", ibeta_);
  vector<double> beta_max;
  parseVector("BETA_MAX", beta_max);
  // check dimension of beta_max
  if(beta_max.size()!=static_cast<unsigned>(getNumberOfArguments()))
    error("Size of BETA_MAX array should be the same as of the number of arguments in ARG");

  // number of bins for beta ladder
  unsigned nbin;
  parse("NBIN", nbin);
  // calculate exponents
  double ibeta_max = static_cast<double>(nbin);
  for(unsigned i=0; i<beta_max.size(); ++i){
    expo_.push_back(std::log(beta_max[i])/std::log(ibeta_max));
  }

  parse("W0", w0_);
  parse("BIASFACTOR", biasf_);
  // allocate beta grid and set bias to zero
  for(unsigned i=0; i<nbin; ++i){
    // bias grid
    bias_.push_back(0.0);
    // beta ladder
    double beta = exp( static_cast<double>(i) / static_cast<double>(nbin-1) * std::log(ibeta_max) );
    beta_.push_back(beta);
  } 
  // print bias to file
  parse("BSTRIDE", Biasstride_);
  parse("BFILE",   Biasfilename_);
  
  // share across replicas or not
  parseVector("SHARED", shared_);

  // monte carlo stuff
  parse("MC_STEPS",MCsteps_);
  parse("MC_STRIDE",MCstride_);
  // adjust for multiple-time steps
  MCstride_ *= getStride();

  // get temperature
  double temp=0.0;
  parse("TEMP",temp);
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();

  checkRead();

  log.printf("  temperature of the system in energy unit %f\n",kbt_);
  log.printf("  initial value of beta %f\n",beta_[ibeta_]);
  log.printf("  number of bins in beta grid %u\n",nbin);
  log.printf("  biasfactor %f\n",biasf_);
  log.printf("  initial hills height %f\n",w0_);
  log.printf("  stride to write bias to file %u\n",Biasstride_);
  log.printf("  write bias to file : %s\n",Biasfilename_.c_str());
  log.printf("  number of replicas %u\n",nrep_);
  log.printf("  number of MC steps %d\n",MCsteps_);
  log.printf("  do MC every %d steps\n", MCstride_);
  log.printf("\n");

  // add components
  addComponent("ibeta");    componentIsNotPeriodic("ibeta");
  addComponent("accbeta");  componentIsNotPeriodic("accbeta");
  addComponent("wtbias");   componentIsNotPeriodic("wtbias");
  
  // initialize random seed
  srand (time(NULL));
  
  // read bias if restarting
  if(getRestart()) read_bias();
}

Rescale::~Rescale()
{
  Biasfile_.close();
}

void Rescale::read_bias()
{
 double MDtime;
 // open file
 IFile *ifile = new IFile();
 ifile->link(*this);
 if(ifile->FileExist(Biasfilename_)){
    ifile->open(Biasfilename_);
    // read all the lines, store last value of bias
    while(ifile->scanField("MD_time",MDtime)){
     for(unsigned i=0; i<bias_.size(); ++i){
      // convert i to string
      stringstream ss;
      ss << i;
      // label
      string label = "b" + ss.str();
      // read entry
      ifile->scanField(label, bias_[i]);
     }
     // new line
     ifile->scanField();
    }
    ifile->close();
 } else {
    error("Cannot find bias file "+Biasfilename_+"\n"); 
 }
 delete ifile;
}

unsigned Rescale::proposeMove(unsigned x, unsigned xmin, unsigned xmax)
{
 int xmin_i = static_cast<int>(xmin);
 int xmax_i = static_cast<int>(xmax);
 int dx;
 int r = rand() % 2;
 if( r % 2 == 0 ) dx = +1;
 else             dx = -1;
 // new index, integer
 int x_new = static_cast<int>(x) + dx;
 // check boundaries
 if(x_new >= xmax_i) x_new = xmax_i-1;
 if(x_new <  xmin_i) x_new = xmin_i;
 return static_cast<unsigned>(x_new);
}

bool Rescale::doAccept(double oldE, double newE, double oldB, double newB)
{
  bool accept = false;
  // calculate delta energy 
  double delta = ( newE + newB - oldE - oldB) / kbt_;
  // if delta is negative always accept move
  if( delta < 0.0 ){ 
   accept = true;
  }else{
   // otherwise extract random number   
   double s = static_cast<double>(rand()) / RAND_MAX;
   if( s < exp(-delta) ) { accept = true; }
  }
  return accept;
}

void Rescale::doMonteCarlo(long int step, vector<double> args)
{
 bool accept;
 double oldE = 0.0;
 // calculate old energy
 for(unsigned j=0; j<args.size(); ++j){
    // calculate energy term
    double fact = 1.0/pow(beta_[ibeta_], expo_[j]) - 1.0;
    oldE += args[j] * fact;
 } 
 // cycle on MC steps 
 for(unsigned i=0;i<MCsteps_;++i){
  // propose move in ibeta
  unsigned new_ibeta = proposeMove(ibeta_, 0, beta_.size());
  // calculate new energy
  double newE = 0.0;
  for(unsigned j=0; j<args.size(); ++j){
    // calculate energy term
    double fact = 1.0/pow(beta_[new_ibeta], expo_[j]) - 1.0;
    newE += args[j] * fact;
  }
  // accept or reject
  accept = doAccept(oldE, newE, bias_[ibeta_], bias_[new_ibeta]);
  if(accept){
   ibeta_ = new_ibeta;
   oldE = newE;
   MCaccbeta_++;
  }
 }
 // send values of parameters to all replicas
 if(comm.Get_rank()==0){
   if(multi_sim_comm.Get_rank()!=0) ibeta_ = 0;
   multi_sim_comm.Sum(&ibeta_, 1); 
 } else {
   ibeta_ = 0;
 }
 // wait for things to be done
 multi_sim_comm.Barrier();
 // local communication
 if(comm.Get_size()>1) comm.Sum(&ibeta_, 1);
 
 // add well-tempered like bias
 double kbDT = kbt_ * ( biasf_ - 1.0 );
 bias_[ibeta_] += w0_ * exp(-bias_[ibeta_] / kbDT);
}

void Rescale::print_bias(long int step)
{
 // if first time open the file
 if(first_bias_){
  first_bias_ = false;
  Biasfile_.link(*this);
  Biasfile_.open(Biasfilename_);
  Biasfile_.setHeavyFlush();
  Biasfile_.fmtField("%30.5f");
 }

 // write fields
 double MDtime = static_cast<double>(step)*getTimeStep();
 Biasfile_.printField("MD_time", MDtime);
 for(unsigned i=0; i<bias_.size(); ++i){
   // convert i to string
   stringstream ss;
   ss << i;
   // label
   string label = "b" + ss.str();
   // print entry
   Biasfile_.printField(label, bias_[i]);
 }
 Biasfile_.printField();
}

void Rescale::calculate()
{
  vector<double> args(getNumberOfArguments(), 0.0);
  // first calculate arguments
  for(unsigned i=0; i<args.size(); ++i){
     args[i] = getArgument(i);
     // sum shared arguments across replicas
     if(shared_[i]==1){
        if(comm.Get_rank()==0) multi_sim_comm.Sum(&args[i], 1);
        else                   args[i] = 0.0;
        if(comm.Get_size()>1)  comm.Sum(&args[i], 1);
     }
  }
    
  // get time step 
  long int step = getStep();
  // do MC at the right time step
  if(step%MCstride_==0&&!getExchangeStep()) doMonteCarlo(step, args);

  // calculate energy and forces
  double ene = 0.0;
  for(unsigned i=0; i<args.size(); ++i){
    // calculate energy term
    double fact = 1.0/pow(beta_[ibeta_], expo_[i]) - 1.0;
    ene += args[i] * fact;
    // add force
    setOutputForce(i, -fact);
  }
  
  // set value of the bias
  setBias(ene);
  // set value of the wt-bias
  getPntrToComponent("wtbias")->set(bias_[ibeta_]);
  // set values of beta
  getPntrToComponent("ibeta")->set(ibeta_);
  // beta acceptance
  if(MCfirst_==-1) MCfirst_=step;
  // calculate acceptance
  double MCtrials = std::floor(static_cast<double>(step-MCfirst_) / static_cast<double>(MCstride_))+1.0;
  double accbeta = static_cast<double>(MCaccbeta_) / static_cast<double>(MCsteps_) / MCtrials;
  getPntrToComponent("accbeta")->set(accbeta);

  // print bias
  if(step%Biasstride_==0) print_bias(step);

}


}
}

