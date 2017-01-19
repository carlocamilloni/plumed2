/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016 The plumed team
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
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Matrix.h"
#include "core/SetupMolInfo.h"
#include "core/ActionSet.h"
#include "tools/File.h"

#include <string>
#include <cmath>
#include <map>
#include <numeric>
#include <ctime>

using namespace std;

namespace PLMD{
namespace isdb{

//+PLUMEDOC COLVAR EM3DSigma
/*
Put Documentation here 

*/
//+ENDPLUMEDOC
   
class EM3DSigma : public Colvar {

private:

 // temperature in kbt
 double kbt_;
 // model GMM - weights and atom types
 vector<double>   GMM_m_w_;
 vector<unsigned> GMM_m_type_;
 // data GMM - means, weights, and covariances
 vector<Vector>             GMM_d_m_;
 vector<double>             GMM_d_w_;
 vector< VectorGeneric<6> > GMM_d_cov_;
 // overlaps 
 vector<double> ovmd_;
 vector<double> ovdd_;
 double ov_cut_;
 vector<double> ovdd_cut_;
 // and derivatives
 vector<Vector> ovmd_der_;
 // constant quantities;
 double cfact_;
 double inv_sqrt2_, sqrt_2pi_;
 // metainference
 unsigned nrep_;
 unsigned replica_;
 vector<double> sigma_mean_;
 // non-marginal version
 vector<double> sigma_;
 // calculate average sigma
 vector<double> ratios_;
 vector<double> Averatios_;
 double nframe_;
 vector<Value*> valueAvesigma_;
 
 // auxiliary stuff
 // list of atom sigmas
 vector<double> s_map_;
 // list of prefactors for overlap between two components of model and data GMM
 // fact_md = w_m * w_d / (2pi)**1.5 / sqrt(det_md)
 vector< double > fact_md_;
 // inverse of the sum of model and data covariances matrices
 vector< VectorGeneric<6> > inv_cov_md_;
 // neighbor list
 double   nl_cutoff_;
 unsigned nl_stride_;
 bool first_time_, no_aver_;
 vector < unsigned > nl_;
 // parallel stuff
 bool serial_;
 unsigned size_;
 unsigned rank_;
 
 // calculate model GMM weights and covariances - these are constants
 void get_GMM_m(vector<AtomNumber> &atoms);
 // read data GMM file
 void get_GMM_d(string gmm_file);
 // normalize GMM
 void normalize_GMM(vector<double> &w);
 // check GMM data
 void check_GMM_d(VectorGeneric<6> &cov, double w);
 
 // get auxiliary stuff
 void get_auxiliary_stuff();
 // get cutoff in overlap
 void get_cutoff_ov();
 // get fact_md and inv_cov_md
 double get_prefactor_inverse (const VectorGeneric<6> &GMM_cov_0, const VectorGeneric<6> &GMM_cov_1,
        double &GMM_w_0, double &GMM_w_1, 
        VectorGeneric<6> &sum, VectorGeneric<6> &inv_sum);
 // calculate self overlaps between data GMM components - ovdd_
 double get_self_overlap(unsigned id);
 // calculate overlap between two components
 double get_overlap(const Vector &m_m, const Vector &d_m, double &fact_md,
                    const VectorGeneric<6> &inv_cov_md, Vector &ov_der);
 double get_overlap(const Vector &m_m, const Vector &d_m, double &fact_md,
                    const VectorGeneric<6> &inv_cov_md);
 // update the neighbor list
 void update_neighbor_list();
 // calculate overlap
 void calculate_overlap();
  
public:
  static void registerKeywords( Keywords& keys );
  explicit EM3DSigma(const ActionOptions&);
// active methods:
  void prepare();
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(EM3DSigma,"EM3DSigma")

void EM3DSigma::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","atoms for which we calculate the density map");
  keys.add("compulsory","GMM_FILE","file with the parameters of the GMM components");
  keys.add("compulsory","TEMP","temperature");
  keys.addFlag("SERIAL",false,"perform the calculation in serial - for debug purpose");
  keys.addFlag("NO_AVER",false,"don't do ensemble averaging");
  keys.addFlag("RELATIVE",false,"calculate relative error");
  keys.add("compulsory","NL_CUTOFF","The cutoff in overlap for the neighbor list");
  keys.add("compulsory","NL_STRIDE","The frequency with which we are updating the neighbor list");
  keys.add("compulsory","SIGMA","values of the sigma parameter");
  keys.add("compulsory","SIGMA_MEAN","starting value for the uncertainty in the mean estimate");
  useCustomisableComponents(keys);
  keys.addOutputComponent("Avesigma", "default","probability ratio for data GMM components");
}

EM3DSigma::EM3DSigma(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
inv_sqrt2_(0.707106781186548),
sqrt_2pi_(2.506628274631001),
nframe_(0.0), nl_cutoff_(-1.0), nl_stride_(0),
first_time_(true), no_aver_(false), serial_(false)
{
  
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  
  string GMM_file;
  parse("GMM_FILE",GMM_file);
   
  // uncertainty stuff
  double sigma_mean; 
  parse("SIGMA_MEAN",sigma_mean);
  parseVector("SIGMA",sigma_);
  
  // temperature
  double temp=0.0;
  parse("TEMP",temp);
  // convert temp to kbt
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();
 
  // neighbor list stuff
  parse("NL_CUTOFF",nl_cutoff_);
  if(nl_cutoff_<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
  parse("NL_STRIDE",nl_stride_);
  if(nl_stride_<=0) error("NL_STRIDE should be explicitly specified and positive");
  
  // serial or parallel
  parseFlag("SERIAL",serial_);
  if(serial_){
    size_=1; rank_=0;
  } else {
    size_=comm.Get_size(); rank_=comm.Get_rank();
  }
 
  parseFlag("NO_AVER",no_aver_);
 
  checkRead();
  
  // get number of replicas
  if(comm.Get_rank()==0) {
    if(no_aver_){
     nrep_ = 1;
     replica_ = 0;
    } else {
     nrep_ = multi_sim_comm.Get_size();
     replica_ = multi_sim_comm.Get_rank();
    }
  } else {
    nrep_ = 0;
    replica_ = 0;
  }
  comm.Sum(&nrep_,1);
  comm.Sum(&replica_,1);
  
  // divide sigma_mean by the square root of the number of replicas
  sigma_mean /= sqrt(static_cast<double>(nrep_));
 
  log.printf("  atoms involved : ");
  for(unsigned i=0;i<atoms.size();++i) log.printf("%d ",atoms[i].serial());
  log.printf("\n");
  log.printf("  GMM data file : %s\n", GMM_file.c_str());
  if(serial_) log.printf("  serial calculation\n");
  if(no_aver_) log.printf("  without ensemble averaging\n");
  log.printf("  neighbor list overlap cutoff : %lf\n", nl_cutoff_);
  log.printf("  neighbor list stride : %u\n",  nl_stride_);
  log.printf("  uncertainty in the mean estimate %f\n",sigma_mean);
  log.printf("  temperature of the system in energy unit %f\n",kbt_);
  log.printf("  number of replicas %u\n",nrep_);
   
  log<<"  Bibliography "<<plumed.cite("Bonomi, Camilloni, Cavalli, Vendruscolo, Sci. Adv. 2, e150117 (2016)");

  // set constant quantity before calculating stuff
  cfact_ = 1.0/pow( 2.0*pi, 1.5 );
  
  // calculate model GMM constant parameters
  get_GMM_m(atoms);

  // read data GMM parameters
  get_GMM_d(GMM_file);
  log.printf("  number of GMM components : %u\n", static_cast<unsigned>(GMM_d_m_.size()));
  
  // normalize GMMs
  normalize_GMM(GMM_m_w_);
  normalize_GMM(GMM_d_w_);
 
  // get self overlaps between data GMM components
  for(unsigned i=0;i<GMM_d_w_.size();++i) {
      double ov = get_self_overlap(i);
      ovdd_.push_back(ov);
      sigma_mean_.push_back(sigma_mean*ov);
  }
 
  // calculate auxiliary stuff
  get_auxiliary_stuff();
  
  // get cutoff for overlap calculation - avoid millions of exp calculations
  get_cutoff_ov();
  
  // and prepare temporary vectors
  ovmd_.resize(GMM_d_w_.size());
  ratios_.resize(ovmd_.size()*sigma_.size(), 0.0);
  Averatios_.resize(ovmd_.size()*sigma_.size(), 0.0);

  // clear things that are not needed anymore
  GMM_d_cov_.clear();

  // request the atoms
  requestAtoms(atoms);
  
  // prepare components for output
  for(unsigned i=0; i<GMM_d_w_.size(); ++i){
     std::string num; Tools::convert(i,num);
     addComponent("Avesigma_"+num); componentIsNotPeriodic("Avesigma_"+num);
     valueAvesigma_.push_back(getPntrToComponent("Avesigma_"+num));
  }
  
  
}

void EM3DSigma::get_GMM_m(vector<AtomNumber> &atoms)
{
  vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();

  // map of atom types to A and B coefficients of scattering factor
  // f(s) = A * exp(-B*s**2)
  // B is in Angstrom squared
  map<string, double> w_map;
  w_map["C"] = 2.49982; // type 0
  w_map["O"] = 1.97692;  // type 1
  w_map["N"] = 2.20402; // type 2
  w_map["S"] = 5.14099;  // type 3
  // map between an atom type and an index
  map<string, unsigned> type_map;
  type_map["C"]=0;
  type_map["O"]=1;
  type_map["N"]=2;
  type_map["S"]=3;
  // fill in the sigma vector
  s_map_.push_back(15.146);
  s_map_.push_back(8.59722);
  s_map_.push_back(11.1116);
  s_map_.push_back(15.8952);
  
  // check if MOLINFO line is present 
  if( moldat.size()==1 ){
    log<<"  MOLINFO DATA found, using proper atom names\n";
    for(unsigned i=0;i<atoms.size();++i){
      // get atom name
      string name = moldat[0]->getAtomName(atoms[i]);
      char type;
      // get atom type
      char first = name.at(0);
      // GOLDEN RULE: type is first letter, if not a number
      if (!isdigit(first)){
         type = first;
      // otherwise is the second
      } else {
         type = name.at(1);
      }
      // check if key in map
      std::string type_s = std::string(1,type);
      if(w_map.find(type_s) != w_map.end()){
        // save atom type
        GMM_m_type_.push_back(type_map[type_s]);
        // this will be normalized to 1 in the final density
        GMM_m_w_.push_back(w_map[type_s]); 
      } else {
        error("Wrong atom type "+type_s+" from atom name "+name+"\n"); 
      }
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
}

	
void EM3DSigma::check_GMM_d(VectorGeneric<6> &cov, double w)
{
 
 // check if positive defined, by calculating the 3 leading principal minors
 double pm1 = cov[0]; 
 double pm2 = cov[0]*cov[3]-cov[1]*cov[1];
 double pm3 = cov[0]*(cov[3]*cov[5]-cov[4]*cov[4])-cov[1]*(cov[1]*cov[5]-cov[4]*cov[2])+cov[2]*(cov[1]*cov[4]-cov[3]*cov[2]);
 // apply Sylvester’s criterion
 if(pm1<=0.0 || pm2<=0.0 || pm3<=0.0)
  error("check data GMM: covariance matrix is not positive defined");
 
 // check if weight is positive
 if(w<0.0) error("check data GMM: weight must be positive");
}


// read GMM data file in PLUMED format:
void EM3DSigma::get_GMM_d(string GMM_file)
{
 int idcomp, beta;
 double w, m0, m1, m2;
 VectorGeneric<6> cov;
 // open file
 IFile *ifile = new IFile();
 if(ifile->FileExist(GMM_file)){
    ifile->open(GMM_file);
    while(ifile->scanField("Id",idcomp)){
     ifile->scanField("Weight",w);
     ifile->scanField("Mean_0",m0);
     ifile->scanField("Mean_1",m1);
     ifile->scanField("Mean_2",m2);
     ifile->scanField("Cov_00",cov[0]);
     ifile->scanField("Cov_01",cov[1]);
     ifile->scanField("Cov_02",cov[2]);
     ifile->scanField("Cov_11",cov[3]);
     ifile->scanField("Cov_12",cov[4]);
     ifile->scanField("Cov_22",cov[5]);
     ifile->scanField("Beta",beta);
     // check input
     check_GMM_d(cov, w);
     // center of the Gaussian
     GMM_d_m_.push_back(Vector(m0,m1,m2));
     // covariance matrix
     GMM_d_cov_.push_back(cov);
     // weights
     GMM_d_w_.push_back(w);
     // new line
     ifile->scanField();
    }
    ifile->close();
 } else {
    error("Cannot find GMM_FILE "+GMM_file+"\n"); 
 }
 delete ifile;

}

// normalize GMM to sum to 1
// since all the GMM components are individually normalized, we just need to 
// divide each weight for the sum of the weights
void EM3DSigma::normalize_GMM(vector<double> &w)
 {
   double norm = accumulate(w.begin(), w.end(), 0.0);
   for(unsigned i=0; i<w.size(); ++i) w[i] /= norm;
 }
 
 void EM3DSigma::get_auxiliary_stuff()
 {
  VectorGeneric<6> cov, sum, inv_sum;
  // cycle on all atoms types
  for(unsigned i=0; i<4; ++i){
   // the Gaussian in density (real) space is the FT of scattering factor
   // f(r) = A * (pi/B)**1.5 * exp(-pi**2/B*r**2)
   double s = sqrt ( 0.5 * s_map_[i] ) / pi * 0.1;
   // covariance matrix for spherical Gaussian
   cov[0]=s*s; cov[1]=0.0; cov[2]=0.0;
               cov[3]=s*s; cov[4]=0.0;
                           cov[5]=s*s;
   // cycle on all data GMM
   for(unsigned j=0; j<GMM_d_m_.size(); ++j){
     // we need the sum of the covariance matrices
     for(unsigned k=0; k<6; ++k) sum[k] = cov[k] + GMM_d_cov_[j][k];
     // and to calculate its determinant
     double det = sum[0]*(sum[3]*sum[5]-sum[4]*sum[4]);
           det -= sum[1]*(sum[1]*sum[5]-sum[4]*sum[2]);
           det += sum[2]*(sum[1]*sum[4]-sum[3]*sum[2]);
     // the constant part of the prefactor is
     double pre_fact =  cfact_ / sqrt(det);
     // and its inverse
     inv_sum[0] = (sum[3]*sum[5] - sum[4]*sum[4])/det;
     inv_sum[1] = (sum[2]*sum[4] - sum[1]*sum[5])/det;
     inv_sum[2] = (sum[1]*sum[4] - sum[2]*sum[3])/det;
     inv_sum[3] = (sum[0]*sum[5] - sum[2]*sum[2])/det;
     inv_sum[4] = (sum[2]*sum[1] - sum[0]*sum[4])/det;
     inv_sum[5] = (sum[0]*sum[3] - sum[1]*sum[1])/det;
     // now we store the pre_fact
     fact_md_.push_back(pre_fact);
     // and the inverse of the sum
     inv_cov_md_.push_back(inv_sum);    
   } 
  } 
 
 }

// get prefactors
double EM3DSigma::get_prefactor_inverse
(const VectorGeneric<6> &GMM_cov_0, const VectorGeneric<6> &GMM_cov_1,
 double &GMM_w_0, double &GMM_w_1, 
 VectorGeneric<6> &sum, VectorGeneric<6> &inv_sum)
{
 // we need the sum of the covariance matrices
 for(unsigned k=0; k<6; ++k) sum[k] = GMM_cov_0[k] + GMM_cov_1[k];
  
 // and to calculate its determinant
 double det = sum[0]*(sum[3]*sum[5]-sum[4]*sum[4]);
       det -= sum[1]*(sum[1]*sum[5]-sum[4]*sum[2]);
       det += sum[2]*(sum[1]*sum[4]-sum[3]*sum[2]);
      
 // the prefactor is 
 double pre_fact =  cfact_ / sqrt(det) * GMM_w_0 * GMM_w_1;

 // and its inverse
 inv_sum[0] = (sum[3]*sum[5] - sum[4]*sum[4])/det;
 inv_sum[1] = (sum[2]*sum[4] - sum[1]*sum[5])/det;
 inv_sum[2] = (sum[1]*sum[4] - sum[2]*sum[3])/det;
 inv_sum[3] = (sum[0]*sum[5] - sum[2]*sum[2])/det;
 inv_sum[4] = (sum[2]*sum[1] - sum[0]*sum[4])/det;
 inv_sum[5] = (sum[0]*sum[3] - sum[1]*sum[1])/det;

 // return pre-factor
 return pre_fact;
}

double EM3DSigma::get_self_overlap(unsigned id)
{
 vector<double> ov;
 VectorGeneric<6> sum, inv_sum;
 Vector ov_der;
 // start loop
 for(unsigned i=0; i<GMM_d_w_.size(); ++i){
   // call auxiliary method
   double pre_fact = get_prefactor_inverse(GMM_d_cov_[id], GMM_d_cov_[i], 
                                             GMM_d_w_[id],   GMM_d_w_[i], sum, inv_sum); 
   // calculate overlap
   double ov_tmp = get_overlap(GMM_d_m_[id], GMM_d_m_[i], pre_fact, inv_sum, ov_der);
   // add to list
   ov.push_back(ov_tmp);
 }
 // calculate total
 double ov_tot = accumulate(ov.begin(), ov.end(), 0.0);
 // sort in ascending order
 std::sort(ov.begin(), ov.end());
 // get cutoff = nl_cutoff_ * ov_tot
 double ov_cut = ov_tot * nl_cutoff_;
 // integrate tail of ov
 double ov_sum = 0.0;
 for(unsigned i=1; i<ov.size(); ++i){
    ov_sum += ov[i];
    if(ov_sum >= ov_cut){
       ov_cut = ov[i-1];
       break;
    } 
 }
 // store 
 ovdd_cut_.push_back(ov_cut);
 // and return it
 return ov_tot;
}

// this is to avoid the calculation of millions of exp function
// when updating the neighbor list using calculate_overlap
void EM3DSigma::get_cutoff_ov()
{
  // temporary stuff
  unsigned GMM_d_w_size = GMM_d_w_.size();
  // set ov_cut_ to a huge number
  ov_cut_ = 1.0+9;  
  // calculate minimum value needed for cutoff
  for(unsigned i=0; i<GMM_d_w_.size(); ++i){
   for(unsigned j=0; j<GMM_m_w_.size(); ++j){
     // get atom type
     unsigned jtype = GMM_m_type_[j];
     // get index in auxiliary lists
     unsigned kaux = jtype * GMM_d_w_size + i;
     // get prefactor and multiply by weights
     double pre_fact = fact_md_[kaux] * GMM_d_w_[i] * GMM_m_w_[j];
     // calculate ov
     double ov = ovdd_cut_[i] / pre_fact;
     // check
     if(ov < ov_cut_) ov_cut_ = ov;
   }
  }
  // set cutoff
  ov_cut_ = -2.0 * std::log(ov_cut_);
}

// version with derivatives
double EM3DSigma::get_overlap(const Vector &m_m, const Vector &d_m, double &fact_md,
                            const VectorGeneric<6> &inv_cov_md, Vector &ov_der)
{
  // calculate vector difference m_m-d_m
  double md_x = m_m[0] - d_m[0];
  double md_y = m_m[1] - d_m[1];
  double md_z = m_m[2] - d_m[2];
  // calculate product of transpose of md and inv_cov_md
  double p_x = md_x*inv_cov_md[0]+md_y*inv_cov_md[1]+md_z*inv_cov_md[2];
  double p_y = md_x*inv_cov_md[1]+md_y*inv_cov_md[3]+md_z*inv_cov_md[4];
  double p_z = md_x*inv_cov_md[2]+md_y*inv_cov_md[4]+md_z*inv_cov_md[5];
  // calculate product of prod and md
  double ov = md_x*p_x+md_y*p_y+md_z*p_z; 
  // final calculation
  ov = fact_md * exp(-0.5*ov);
  // derivatives
  double x = md_x*inv_cov_md[0] + md_y*inv_cov_md[1] + md_z*inv_cov_md[2];
  double y = md_x*inv_cov_md[1] + md_y*inv_cov_md[3] + md_z*inv_cov_md[4];
  double z = md_x*inv_cov_md[2] + md_y*inv_cov_md[4] + md_z*inv_cov_md[5];
  ov_der = ov * Vector(x, y, z); 
  return ov;
}

// fast version without derivatives and cutoff used for neighbor list
double EM3DSigma::get_overlap(const Vector &m_m, const Vector &d_m, double &fact_md, 
                            const VectorGeneric<6> &inv_cov_md)
                        
{
  // calculate vector difference m_m-d_m
  double md_x = m_m[0] - d_m[0];
  double md_y = m_m[1] - d_m[1];
  double md_z = m_m[2] - d_m[2];
  // calculate product of transpose of md and inv_cov_md
  double p_x = md_x*inv_cov_md[0]+md_y*inv_cov_md[1]+md_z*inv_cov_md[2];
  double p_y = md_x*inv_cov_md[1]+md_y*inv_cov_md[3]+md_z*inv_cov_md[4];
  double p_z = md_x*inv_cov_md[2]+md_y*inv_cov_md[4]+md_z*inv_cov_md[5];
  // calculate product of prod and md
  double ov = md_x*p_x+md_y*p_y+md_z*p_z; 
  // final calculation
  if( ov > ov_cut_ ){ 
    ov = 0.0;
  } else { 
    ov = fact_md * exp(-0.5*ov);
  }
  return ov;
}

void EM3DSigma::update_neighbor_list()
{
  // temp stuff
  unsigned GMM_d_w_size = GMM_d_w_.size();
  unsigned GMM_m_w_size = GMM_m_w_.size();
  // local neighbor list
  vector < unsigned > nl_l;
  // clear old neighbor list
  nl_.clear();
  // cycle on all overlaps (in parallel)
  unsigned nover = GMM_d_w_size * GMM_m_w_size;
  for(unsigned k=rank_; k<nover; k=k+size_){
      // get indexes
      unsigned i = k / GMM_m_w_size;
      unsigned j = k % GMM_m_w_size;
      // get atom type
      unsigned jtype = GMM_m_type_[j];
      // get index in auxiliary lists
      unsigned kaux = jtype * GMM_d_w_size + i;
      // get prefactor and multiply by weights
      double pre_fact = fact_md_[kaux] * GMM_d_w_[i] * GMM_m_w_[j];
      // calculate overlap
      double ov = get_overlap(GMM_d_m_[i], getPosition(j), pre_fact, inv_cov_md_[kaux]);
      // fill the neighbor list
      if(ov >= ovdd_cut_[i]) nl_l.push_back(k);
  }
  // find total dimension of neighborlist
  vector <int> recvcounts(size_, 0);
  recvcounts[rank_] = nl_l.size();
  comm.Sum(&recvcounts[0], size_);
  int tot_size = accumulate(recvcounts.begin(), recvcounts.end(), 0);
  // resize neighbor stuff
  nl_.resize(tot_size);
  // calculate vector of displacement
  vector<int> disp(size_);
  disp[0] = 0;
  int rank_size = 0;
  for(unsigned i=0; i<size_-1; ++i){
    rank_size += recvcounts[i];
    disp[i+1] = rank_size;
  }
  // Allgather neighbor list
  comm.Allgatherv(&nl_l[0], recvcounts[rank_], &nl_[0], &recvcounts[0], &disp[0]);
  // now resize derivatives
  ovmd_der_.resize(tot_size);
}

void EM3DSigma::prepare()
{
  if(getExchangeStep()) first_time_=true;
}


// overlap calculator
void EM3DSigma::calculate_overlap(){

  //makeWhole();
  if(first_time_ || getExchangeStep() || getStep()%nl_stride_==0){
     update_neighbor_list();
     first_time_=false;
  }
  
  // clean temporary vectors
  for(unsigned i=0; i<ovmd_.size(); ++i)     ovmd_[i] = 0.0;
  for(unsigned i=0; i<ovmd_der_.size(); ++i) ovmd_der_[i] = Vector(0,0,0);
  
  // we have to cycle over all model and data GMM components in the neighbor list
  unsigned GMM_d_w_size = GMM_d_w_.size();
  unsigned GMM_m_w_size = GMM_m_w_.size();
  for(unsigned i=rank_;i<nl_.size();i=i+size_) {
      // get indexes of data and model component
      unsigned id = nl_[i] / GMM_m_w_size;
      unsigned im = nl_[i] % GMM_m_w_size;
      // get atom type
      unsigned jtype = GMM_m_type_[im];
      // get index in auxiliary lists
      unsigned kaux = jtype * GMM_d_w_size + id;
      // get prefactor and multiply by weights
      double pre_fact = fact_md_[kaux] * GMM_d_w_[id] * GMM_m_w_[im];
      // add overlap with im component of model GMM
      ovmd_[id] += get_overlap(GMM_d_m_[id], getPosition(im), pre_fact,
                               inv_cov_md_[kaux], ovmd_der_[i]);
  }
  // if parallel, communicate stuff
  if(!serial_){
   comm.Sum(&ovmd_[0], ovmd_.size());
   comm.Sum(&ovmd_der_[0][0], 3*ovmd_der_.size());
  }
}


void EM3DSigma::calculate(){

  // calculate CV 
  calculate_overlap();
  
  // rescale factor for ensemble average
  double escale = 1.0 / static_cast<double>(nrep_);
  
  // calculate average of ovmd_ across replicas
  if(!no_aver_){
   if(comm.Get_rank()==0){
      multi_sim_comm.Sum(&ovmd_[0], ovmd_.size());
      for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_[i] *= escale;
   } else {
      for(unsigned i=0; i<ovmd_.size(); ++i) ovmd_[i]  = 0.0;
   }
   comm.Sum(&ovmd_[0], ovmd_.size());
  }
 
  // clear ratios_
  for(unsigned i=0; i<ratios_.size(); ++i) ratios_[i] = 0.0;
 
  // calculate marginal and non marginal probabilities and ratio
  for(unsigned k=rank_; k<ratios_.size(); k=k+size_){
     // get indexes
     unsigned i = k / sigma_.size();
     unsigned j = k % sigma_.size();
     // calculate  err function
     double err_f = erf ( ( ovmd_[i]-ovdd_[i] ) * inv_sqrt2_ / sigma_mean_[i] );
     // marginal probability
     double marg_p = 0.5 / (ovmd_[i]-ovdd_[i]) * err_f;
     // calculate (absolute) effective sigma
     double sigma_eff2 =  sigma_mean_[i]*sigma_mean_[i] + sigma_[j]*sigma_[j]*ovdd_[i]*ovdd_[i];
     // non-marginal probability
     double non_marg_p = 1.0 / sqrt_2pi_ / sigma_eff2 * exp( - 0.5 * ( ovmd_[i]-ovdd_[i] ) * ( ovmd_[i]-ovdd_[i] ) / sigma_eff2);
     // calculate ratio of probabilities 
     ratios_[k] = non_marg_p / marg_p; 
  }
  
  // now sum if parallel
  if(!serial_) comm.Sum(&ratios_[0], ratios_.size());
  
  // increment average probabilities over the trajectory
  for(unsigned i=0; i<ratios_.size(); ++i) Averatios_[i] += ratios_[i];
  
  // increment frame counter
  nframe_ += 1.0;
  
  // calculate averages
  vector<double> num(ovmd_.size(), 0.0);
  vector<double> den(ovmd_.size(), 0.0);

  for(unsigned k=rank_; k<ratios_.size(); k=k+size_){
     // get indexes
     unsigned i = k / sigma_.size();
     unsigned j = k % sigma_.size();

     // skip last point on grid
     if(j==sigma_.size()-1) continue;
     
     // calculate p
     double p = Averatios_[k] / nframe_;
     
     // increment num and den
     num[i] += sigma_[j] * p * ( sigma_[j+1] - sigma_[j] );
     den[i] += p * ( sigma_[j+1] - sigma_[j] );  
  }
  
  // now sum if parallel
  if(!serial_){
    comm.Sum(&num[0], num.size());
    comm.Sum(&den[0], den.size());
  }
  
  // set value of components
  for(unsigned i=0; i<valueAvesigma_.size(); ++i) valueAvesigma_[i]->set(num[i]/den[i]);
}

}
}
