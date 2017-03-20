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
 This class was originally written by Alexander Jussupow and
 Carlo Camilloni. Hierarchical algorithm and the middleman algorithm by Max Muehlbauer 
*/
#include <bitset>
#include "Colvar.h"
#include "ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "core/SetupMolInfo.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <string>
#include <cmath>
#include <map>
#include<time.h>
using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR HSAXS
/*

Calculate SAXS scattered intensity

*/
//+ENDPLUMEDOC
   
class HSAXS : public Colvar {
private:
  bool                     pbc;
  bool                     serial;
  unsigned                 numq;
  int      		   truncation;
  vector<double>           q_list;
  vector<vector<double> >  FF_value;
  vector<double>           FF_rank;
  vector<double>           av;
  vector<double>           bv;
  vector<double>           dv;
  vector<double>           rcoeff;
  vector<double>           tcoeff;
  unsigned int             p2;
  double                   maxdist;
  double                   maxq;
  unsigned                 numl;

public:
  static void              registerKeywords( Keywords& keys );
  explicit                 HSAXS(const ActionOptions&);
  void                     getMartiniSFparam(const vector<AtomNumber> &atoms, vector<vector<long double> > &parameter);
  void                     calculateASF(const vector<AtomNumber> &atoms, vector<vector<long double> > &FF_tmp);
  void                     derivatives();
  virtual void             calculate();
  Vector2d                 dXHarmonics(unsigned k, unsigned int i, int n, int m, vector<Vector2d> &decRnm);
  Vector2d                 dYHarmonics(unsigned k, unsigned int i, int n, int m, vector<Vector2d> &decRnm);
  Vector2d                 dZHarmonics(unsigned k, unsigned int i, int n, int m, vector<Vector2d> &decRnm);
  void		           cal_coeff();
  void                     cal_rotcoeff(); 
  void                     cal_transcoeff(int level, double qval, int tr);
  std::bitset<24>          cal_bin(double dec);
  std::bitset<72>          interleave(std::bitset<24> x, std::bitset<24> y, std::bitset<24> z);
  void                     deinterleave(std::bitset<72> index, std::bitset<24> &x, std::bitset<24> &y, std::bitset<24> &z);
  long unsigned            l_index(unsigned l, std::bitset<72> coord);
  Vector                   cal_boxcentre(unsigned index, unsigned l);
  long unsigned            own_index(unsigned ownLevel, double x, double y, double z);
  Vector                   own_centre(unsigned ownLevel, double x, double y, double z);
  long unsigned            parent_index(long unsigned index);
  long unsigned            parent_index(unsigned ownLevel, double x, double y, double z);
  Vector                   parent_centre(unsigned ownLevel, double x, double y, double z);
  double                   factorial(double n, double m); 
};

PLUMED_REGISTER_ACTION(HSAXS,"HSAXS")

void HSAXS::registerKeywords(Keywords& keys){
  Colvar::registerKeywords( keys );
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("ATOMISTIC",false,"calculate SAXS for an atomistic model");
  keys.add("atoms","ATOMS","The atoms to be included in the calculation, e.g. the whole protein.");
  keys.add("compulsory","NUMQ","Number of used q values");
  keys.add("compulsory","NUML","0","Number of levels in the hierarchical octree structure");
  keys.add("numbered","QVALUE","Selected scattering lenghts in Angstrom are given as QVALUE1, QVALUE2, ... .");
  keys.add("numbered","PARAMETERS","Used parameter Keywords like PARAMETERS1, PARAMETERS2. These are used to calculate the structure factor for the i-th atom/bead.");
  keys.addFlag("ADDEXPVALUES",false,"Set to TRUE if you want to have fixed components with the experimental values.");
  keys.add("numbered","EXPINT","Add an experimental value for each q value.");
  keys.add("compulsory","SCEXP","1.0","SCALING value of the experimental data. Usefull to simplify the comparison.");
}
//constructor
HSAXS::HSAXS(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false)
{
  //read in atoms 
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  const unsigned size = atoms.size();

  parseFlag("SERIAL",serial);
  //no pbcs used
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  //read number of qvalues
  numq = 0;
  parse("NUMQ",numq);
  if(numq==0) error("NUMQ must be set");  
  //read in a truncation value
  truncation= 0;
  //set number of levels
  numl=0;
  parse("NUML",numl);
  if(numl==0) numl=0;
  //read experimental scaling factor
  double scexp = 0;
  parse("SCEXP",scexp);
  if(scexp==0) scexp=1.0;  
  //store qvalues in q_list
  q_list.resize( numq );
  unsigned ntarget=0;
  for(unsigned i=0;i<numq;++i){
    if( !parseNumbered( "QVALUE", i+1, q_list[i]) ) break;
    ntarget++;
  }
  if( ntarget!=numq ) error("found wrong number of qvalue values");

  for(unsigned i=0;i<numq;i++) {
    if(q_list[i]==0.) error("it is not possible to set q=0\n"); 
    log.printf("  my q: %lf \n",q_list[i]);    
    //Tools::convert converts qvalues of type double to strings
    std::string num; Tools::convert(i,num);
    addComponentWithDerivatives("q_"+num);
    componentIsNotPeriodic("q_"+num);
  }
  //this flag should switch between MARTINI and atomistic calculations
  bool atomistic=false;
  parseFlag("ATOMISTIC",atomistic);
  vector<vector<long double> >  FF_tmp;
  FF_tmp.resize(numq,vector<long double>(size));
  if(!atomistic) {
    //read in parameter vector for MARTINI parameters
    vector<vector<long double> > parameter;
    parameter.resize(size);
    ntarget=0;
    for(unsigned i=0;i<size;++i){
      if( !parseNumberedVector( "PARAMETERS", i+1, parameter[i]) ) break;
      ntarget++;
    }
    if( ntarget==0 ) getMartiniSFparam(atoms, parameter);
    else if( ntarget!=size ) error("found wrong number of parameter vectors");
    for(unsigned i=0;i<size;++i) {
      for(unsigned k=0;k<numq;++k){
        for(unsigned j=0;j<parameter[i].size();++j) {
          FF_tmp[k][i]+= parameter[i][j]*powl(static_cast<long double>(q_list[k]),j);
        }
      }
    }
  } else {
    //form factors for atomistic calculations
    calculateASF(atoms, FF_tmp);
  }
  // Calculate Rank of FF_matrix and set the entries of the formfactor matrix where k is the qvalue and i the atomnumber
  FF_rank.resize(numq);
  FF_value.resize(numq,vector<double>(size));
  for(unsigned k=0;k<numq;++k){
    for(unsigned i=0;i<size;i++){
       FF_value[k][i] = static_cast<double>(FF_tmp[k][i]);
       FF_rank[k]+=FF_value[k][i]*FF_value[k][i];
    }
  }
  bool exp=false;
  parseFlag("ADDEXPVALUES",exp);
  if(exp) {
    vector<double>   expint;
    expint.resize( numq ); 
    unsigned ntarget=0;
    for(unsigned i=0;i<numq;++i){
       if( !parseNumbered( "EXPINT", i+1, expint[i] ) ) break;
       ntarget++; 
    }
    if( ntarget!=numq ) error("found wrong number of EXPINT values");
    for(unsigned i=0;i<numq;i++) {
      std::string num; Tools::convert(i,num);
      addComponent("exp_"+num);
      componentIsNotPeriodic("exp_"+num);
      Value* comp=getPntrToComponent("exp_"+num); comp->set(expint[i]*scexp);
    }
  }
  // convert units to nm^-1
  maxq=0;
  for(unsigned i=0;i<numq;++i){
    q_list[i]=q_list[i]*10.0;    //factor 10 to convert from A^-1 to nm^-1
    if(q_list[i]>maxq) maxq=q_list[i];
  } 
  requestAtoms(atoms);
  checkRead();
  //precalculates coefficients for derivatives
  cal_coeff();
}


//approximation of the Debye function for SAXS profiling by a truncated harmonic expansion using clever hierarchical data structures (see Gumerov et al. 2012 and Berlin et al. 2014)
void HSAXS::calculate(){
  if(pbc) makeWhole();
  //get parameters for mpi parallelization
  unsigned           stride                 = comm.Get_size();
  unsigned           rank                   = comm.Get_rank();
  if(serial){
    stride                                  = 1;
    rank                                    = 0;
  }
  //number of atoms
  const unsigned int size                   = getNumberOfAtoms();
  //scattering profile I
  vector<double>     I(numq,0);
  //calculate the distance of the centre of the box at lmax for each atom
  vector<Vector>     centreDist(size);
  //coordinates in the unit box. Useful for calculating pertaining box, centre and parent. Compare Gumeroc, 2005, Chapter 5
  vector<Vector>     unitCoords(size);
  Vector             min                    = getPosition(0);
  Vector             max                    = getPosition(0);
  for(unsigned i=rank;i<size;i+=stride) {
   centreDist[i]                            = getPosition(i);
   if(centreDist[i][0]<min[0]) min[0]       = centreDist[i][0];
   if(centreDist[i][1]<min[1]) min[1]       = centreDist[i][1];
   if(centreDist[i][2]<min[2]) min[2]       = centreDist[i][2];
   if(centreDist[i][0]>max[0]) max[0]       = centreDist[i][0];
   if(centreDist[i][1]>max[1]) max[1]       = centreDist[i][1];
   if(centreDist[i][2]>max[2]) max[2]       = centreDist[i][2];
  }
  max                                      -= min;
  maxdist                                   = max[0];
  if(maxdist<max[1]) maxdist                = max[1];
  if(maxdist<max[2]) maxdist                = max[2];
  log<<"min "<<min<<" maxdist "<<maxdist<<"\n";
  for(unsigned i=rank;i<size;i+=stride) {
    centreDist[i]                           = (centreDist[i] - min)/maxdist;
  }
  vector<Vector> polar(size);
  for(unsigned i=rank;i<size;i+=stride) {
    //the coordinates of the atoms in the unit box
    unitCoords[i]                           = centreDist[i];
    centreDist[i]                          -= own_centre(numl,unitCoords[i][0],unitCoords[i][1],unitCoords[i][2]);
    centreDist[i]                          *= maxdist;
    //r
    polar[i][0]                             = sqrt(centreDist[i][0]*centreDist[i][0]+centreDist[i][1]*centreDist[i][1]+centreDist[i][2]*centreDist[i][2]);
    //cos(theta)
    polar[i][1]                             = centreDist[i][2]/polar[i][0];
    //phi
    polar[i][2]                             = atan2(centreDist[i][1],centreDist[i][0]);

  }
  truncation                                = int(maxdist*maxq);
  if(truncation<pow(2.,numl+3.)) truncation = pow(2.,numl+3.);
  //number of terms in the basis set expansion
  p2                                        = truncation*truncation;
  log<<"truncation "<<truncation<<"\n";
  //re-expansion coefficients for translating expansion coefficients by the RCR-formalism. eq. 41, Gumerov 2012
  cal_rotcoeff();
  //as the legndre polynomials and the exponential term in the basis set expansion are not function of the scattering wavenumber q, they can be precomputed here
  vector<Vector2d> qRnm(p2*size);
  for(unsigned int i=rank;i<size;i+=stride) {
    for(int n=0;n<truncation;n+=1) {
      for(int m=0;m<(n+1);m++) {
        int order             = m  - n;
        int x                 = p2 * i + n  * n + m;
        double gsl            = gsl_sf_legendre_sphPlm(n,abs(order),polar[i][1]);
        //real part of the spherical basis function of order m, degree n of atom i
        qRnm[x][0]            = gsl * cos(order*polar[i][2]);
        //imaginary part of the spherical basis function of order m, degree n of atom i
        qRnm[x][1]            = gsl * sin(order*polar[i][2]);
        }
      }
    }
  //sum over qvalues
  for (int k=numq-1; k>=0; k--) {
    int    trunc                = truncation/pow(2.,numl+0.);
    double p22                  = trunc*trunc;
    int    num_boxes            = pow(8.,numl+0.);

    //double sum over the p^2 expansion terms
    vector<Vector2d> Bnm(p22*num_boxes);
    vector<Vector2d> bu         = qRnm;
    for(unsigned int i=rank;i<size;i+=stride) {
      double pq                 = polar[i][0]* q_list[k];
      int index=own_index(numl,unitCoords[i][0],unitCoords[i][1],unitCoords[i][2]);
      for(int n=trunc-1;n>=0;n-=1) {
        //the spherical bessel functions do not depend on the order and are therefore precomputed here
        double bessel           = gsl_sf_bessel_jl(n,pq);
        //here conj(R(m,n))=R(-m,n) is used to decrease the terms in the sum over m by a factor of two
        for(int m=0;m<(n+1);m++) {
          int order             = m-n;
          int s                 = p22*index + n  * n + m;
          int t                 = s  - 2 * order;
          int x                 = p2 * i +  n  * n + m;
          int y                 = p2 * i +  n  * n + m - 2 * order;
          //real part of the spherical basis function of order m, degree n of atom i
          qRnm[x]              *= bessel;
          //real part of the spherical basis function of order -m, degree n of atom i
          qRnm[y][0]            = qRnm[x][0];
          //imaginary part of the spherical basis function of order -m, degree n of atom i
          qRnm[y][1]            =-qRnm[x][1];
          //expansion coefficient of order m and degree n
          Bnm[s]               += FF_value[k][i] * qRnm[y];
          //correction for expansion coefficient of order -m and degree n
          if(order!=0) Bnm[t]  += FF_value[k][i] * qRnm[x];
        }   
      }
    }
  if(!serial) {
    comm.Sum(&Bnm[0][0],2*p22*num_boxes);
  }
  //here the coefficients are caculated for the lowest level; therefore, the upward pass should follow now
  vector<Vector2d>             current                    = Bnm;
  //loop over levels; omp?
    for(unsigned int l=numl;l>=1;l--) {
      int                          old_num_boxes          = num_boxes;
      int                          old_trunc              = trunc;
      int                          old_p22                = p22;
      //becaude maxdist in parent is twice as high and trunc = q *maxdist (+ error_bound)
      trunc                                               = old_trunc*2;
      p22                                                 = trunc*trunc;
      num_boxes                                           = old_num_boxes/8;
      //coefficcients for the translation part of the rcr decomposition. Depend on q and t and therefore cannot be precomputed.
      cal_transcoeff(l,q_list[k],trunc+2);
      vector<Vector2d>             previous               = current;
      current.clear();
      current.resize(p22*num_boxes);
      //loop over boxes: should be mpi parallelized at some point
      for(int i=0;i<old_num_boxes;i++) {
        vector<Vector2d>           boxvector(old_p22);
        for(int j=0;j<old_p22;j++) boxvector[j]           = previous[old_p22*i + j]; 
        //translation vector between the own box centre and the box centre of the parent box
        Vector                     t                      = cal_boxcentre(parent_index(i),l-1)-cal_boxcentre(i,l);
        //euler angles for the rotations in the rcr decomposition of the translation; note, that these are different than the ones mentioned in the appendix of Gumerov et al. Comp. Chem. (2012)
        //especially note, that the sign of alpha is inverted to the one described in eq. A22, Gumerov 2012. No mathematical reason.Just makes the algorithm work for half-boxes.
        //beta is considered indirectly later on
        //gamma is arbitrary, therefore, set to 0 and therefore does not appear at all in this formulation
        double                     alpha                  = -atan2(t[1],t[0]);
        for(int n=0;n<trunc;n++) {
          int                      n4                     = (4*pow(n+0.,3.)-n)/3;
          int                      n2                     = 2*n+1;
          for(int m=0;m<2*n+1;m++) {
            int                    order                  = m-n;
            double                 ap                     = -order*alpha;
            Vector2d               backrot;
            //backrotation. Compare eq. A22 in Gumerov 2012 with exchanged alpha and gamma
            for(int mp=0;mp<2*n+1;mp++) {
              Vector2d             trans;
              int                  orderp                 = mp-n;
              //coaxial translation. Compare eq. A 15 in Gumerov 2012
              for(int np=abs(orderp);np<old_trunc;np++) {
                int                np4                    = (4*pow(np+0.,3.)-np)/3;
                int                np2                    = 2*np+1;
                double             Rmnnp                  = tcoeff[n*(n+1)*(n+2)/6 + abs(mp-n)*(n+1-0.5*(abs(mp-n)-1)) + np-abs(mp-n)];
                if(np>n)           Rmnnp                  = pow(-1.,n+np+0.) * tcoeff[np*(np+1)*(np+2)/6   +   abs(orderp)*(np+1-0.5*(abs(orderp)-1))   +   n-abs(orderp)]; 
                Vector2d           rot;
                //rotation. Compare eq. A22 in Gumerov, 2012
                for(int mpp=0;mpp<2*np+1;mpp++) {
                  double           app                    = (mpp-np)*alpha;
                  //the if clause is in principal an indirect consideration of beta, which makes it unnecessary to calculate a normalization
                  double           Hmmn                   = rcoeff[np4 + np2*(mp-n+np) + mpp];
                  if(t[2]<0.)      Hmmn                   = rcoeff[np4 + np2*(mp-n+np) + mpp-2*(mpp-np)] * pow(-1.,np+orderp+mpp-np+0.);
                  rot[0]                                 += (cos(app) * boxvector[np*np + mpp][0] - sin(app) * boxvector[np*np + mpp][1]) * Hmmn;
                  rot[1]                                 += (cos(app) * boxvector[np*np + mpp][1] + sin(app) * boxvector[np*np + mpp][0]) * Hmmn;
                }
                trans                                    += Rmnnp * rot;
              }
              //the if clause is in principal an indirect consideration of beta, which makes it unnecessary to calculate a normalization
              double               Hmmnb                  = rcoeff[n4 + n2*m + mp];
              if(t[2]<0.)          Hmmnb                  = rcoeff[n4 + n2*m + mp-2*orderp] * pow(-1.,n+order+orderp+0.);
              backrot				         += trans*Hmmnb;
            }
            //the sine and cosine terms come from the exponential term in front of the seum in eq. A22, Gumerov, 2012
            current[p22*parent_index(i) + n  * n + m][0] += (cos(ap) * backrot[0] - sin(ap) * backrot[1]);
            current[p22*parent_index(i) + n  * n + m][1] += (cos(ap) * backrot[1] + sin(ap) * backrot[0]);
          }
        }
      } 
    }
    for(int t=0;t<num_boxes;t++) {
      for(int n=rank;n<trunc;n+=stride) {
        for(int m=0;m<(2*n)+1;m++) {
          int s          = n * n + m+p22*t;
          I[k]          += current[s][0] * current[s][0] + current[s][1] * current[s][1];
        }
      }
    }
    //reset the spherical harmonics to the precalculated legendre part
    qRnm=bu;
  } //end of the k loop
  if(!serial) {
    comm.Sum(&I[0],numq);
  }
  //output scattering profile
  for(unsigned k=0; k<numq; k++) {
    Value* val            = getPntrToComponent(k);
    I[k]                  = 4 * M_PI*I[k];
    log<<"I["<<k<<"] = "<<I[k]<<"\n";
    val->set(I[k]);
  }
}



//calculation of the factorial dependent normalization factor using gamma functions to avoid overflow for large n.
//needed for the calculations of rotational expansion coefficients according to eq. A25 in Gumerov, 2012
double HSAXS::factorial(double n, double m) {
  return tgamma(2.*n+1.)/(tgamma(-m+n+1.)*tgamma(m+n+1.));
}

//coefficients for partial derivatives of the spherical basis functions
void HSAXS::cal_coeff() {
  //adjust the maximum truncation for supplemental coefficients if necessary. However, if you ever need a higher value than 1000, you're probably doing something unsound
  int                               td         = 1000;
  int                               p2d        = td * td;
  av.resize(p2d);
  bv.resize(p2d);
  dv.resize(p2d);
  for( int n=0;n<td;n++) {
    for( int m=0;m<(2*n)+1;m++) {
      double                        mval       = m - n;
      double                        nval       = n;
      av[n*n+m]                                = sqrt(((nval+mval+1)*(nval+1-mval))/(((2*nval)+1)*((2*nval)+3)));
      bv[n*n+m]                                = sqrt(((nval-mval-1)*(nval-mval))/(((2*nval)-1)*((2*nval)+1)));
      if((-n<=(m-n)) && ((m-n)<0))  bv[n*n+m] *= -1;
      dv[n*n+m]                                = 0.5 * sqrt((nval-mval)*(nval+mval+1));
      if(m-n<0)                     dv[n*n+m] *= -1;
    }
  }
}


void HSAXS::cal_rotcoeff() {
  //number of rotational reexpansion coefficients. Evaluate the sum from n=0 to p over (2n+1)^2
  double numrc=(truncation+0.)/3*(4*pow((truncation+0.),2.)-1);
  vector<double> rv;
  //vector of coefficients for beta=pi/2 --> are used to construct the coefficients for arbitrary angles
  rv.resize(numrc);
  //vector of coefficients for arbitrary angles, i.e. beta=acos(1/sqrt(3))
  rcoeff.resize(numrc);
  for(int n=0;n<truncation;n++) {
    //m is used as m'+n or m+n as necessary
    for( int m=n;m<(2*n)+1;m++) {
      //as usual: m is the access index, order is the number used in the math
      int order = m-n;
      //offset for n. Corresponds to the number of elements in an expansion with truncation n.
      int n4    = (4*pow(n+0.,3.)-n)/3;
      //offset for m. Corresponds to the number of m' per m
      int n2    = 2*n+1;
      //coefficients for m'=0 and m=0,... + symmetries (mm'=m'm=-m-m'=-m'-m). See eq. A25 for analytical formula and A27 for symmetries, Gumerov 2012
      //n 0 mp
      rv[n4+n2*n+m]                           = pow(-1.,m-n+0.) * gsl_sf_legendre_sphPlm(n,abs(order),0)/sqrt((2*n+1)/(4*M_PI));
      //n mp 0
      rv[n4+n2*m+n]                           = rv[n4+n2*n+m];
      //n -mp 0
      rv[n4+n2*(m-2*(m-n))+n]                 = rv[n4+n2*n+m];
      //n 0 -mp
      rv[n4+n2*n+m-2*(m-n)]                   = rv[n4+n2*n+m];
      //coefficients for m=n and m'=0,... + symmetries. Analytical formula as given in eq. A25, Gumerov 2012
      //n n mp
      rv[n4+n2*2*n+m]                         = sqrt(factorial(n,order)) * pow(cos(M_PI/4),n+order+0.) * pow(sin(M_PI/4),n-order-0.);
      if(-(m-n)<0) rv[n4+n2*2*n+m]           *= pow(-1.,m-n+0.);
      //n mp n
      rv[n4+n2*m+2*n]                         = rv[n4+n2*2*n+m];
      //n -n -mp
      rv[n4+m-2*(m-n)]                        = rv[n4+n2*2*n+m];
      //n -mp -n
      rv[n4+n2*(m-2*(m-n))]                   = rv[n4+n2*2*n+m];
      //n n -mp
      rv[n4+n2*2*n+m-2*(m-n)]                 = pow(-1.,2.0*n+order) * rv[n4+n2*2*n+m];
      //n -mp n
      rv[n4+n2*(m-2*(m-n))+2*n]               = pow(-1.,2.0*n+order) * rv[n4+n2*2*n+m];
      //n -n mp
      rv[n4+m]                                = pow(-1.,2.0*n+order) * rv[n4+n2*2*n+m];
      //n mp -n
      rv[n4+n2*m]                             = pow(-1.,2.0*n+order) * rv[n4+n2*2*n+m];
    }
  }
  //recursive calculation of  mp=1 from mp=0. eq. A26 with m=0 and beta=pi/2 in Gumerov 2012
  for(int n=2;n<truncation;n++) {
    for( int m=n;m<(2*n);m++) {
      int order = m-n;
      int n4    = (4*pow(n+0.,3.)-n)/3;
      int n2    = 2*n+1;
      int n4p   = (4*pow(n-1.,3.)-(n-1))/3;
      int n2p   = (2*(n-1)+1);
      //recursive calculation of the coefficients for m'=1 and m=0,... from the coefficients of m'=0 and m=0,...
      //n 1 mp
      rv[n4p+n2p*n+m-1]                       = 1/bv[n*n+n] * 0.5 * (bv[n*n+m-1-2*order] * rv[n4+n2*n+m+1] - bv[n*n+m-1] * rv[n4+n2*n+m-1])- 1/bv[n*n+n] * av[(n-1)*(n-1)+m-1] * rv[n4+n2*n+m];
      //n mp 1
      rv[n4p+n2p*(m-1)+n]                     = rv[n4p+n2p*n+m-1];
      //n -1 -mp
      rv[n4p+n2p*(n-2)+m-1-2*(m-n)]           = rv[n4p+n2p*n+m-1];
      //n -mp -1
      rv[n4p+n2p*(m-1-2*(m-n))+n-2]           = rv[n4p+n2p*n+m-1];
      //n -1 mp
      rv[n4p+n2p*(n-2)+m-1]                   = pow(-1.,n-1+m-n+1.) * rv[n4p+n2p*n+m-1];
      //n 1 -mp
      rv[n4p+n2p*n+m-1-2*(m-n)]               = pow(-1.,n-1+m-n+1.) * rv[n4p+n2p*n+m-1]; 
      //n mp -1
      rv[n4p+n2p*(m-1)+n-2]                   = pow(-1.,n-1+m-n+1.) * rv[n4p+n2p*n+m-1];
      //n -mp 1
      rv[n4p+n2p*(m-1-2*(m-n))+n]             = pow(-1.,n-1+m-n+1.) * rv[n4p+n2p*n+m-1];
    }
  }
  //second recursion to fill the matrix. eq. A23 in Gumerov 2012
  for(int n=2;n<truncation;n++) {
    for(int mp=n+1;mp<2*n-1;mp++) {
      for( int m=n+2;m<2*n;m++) {
        //necessary if clause because of compulsory order of for loops
        if(mp<m) {
          int n4 = (4*pow(n+0.,3.)-n)/3;
          int n2 = 2*n+1;
          //n m mp+1
          rv[n4+n2*m+mp+1]                      = -1/dv[n*n+mp] * (dv[n*n+m-1]*rv[n4+n2*(m-1)+mp] - dv[n*n+m]*rv[n4+n2*(m+1)+mp] - dv[n*n+mp-1]*rv[n4+n2*m+mp-1]);
          //n mp+1 m
          rv[n4+n2*(mp+1)+m]                    = rv[n4+n2*m+mp+1];
          //n -(mp+1) -m
          rv[n4+n2*(mp+1-2*(mp+1-n))+m-2*(m-n)] = rv[n4+n2*m+mp+1];
          //n -m -(mp+1)
          rv[n4+n2*(m-2*(m-n))+mp+1-2*(mp+1-n)] = rv[n4+n2*m+mp+1];
          //n -m mp+1
          rv[n4+n2*(m-2*(m-n))+mp+1]            = rv[n4+n2*m+mp+1] *  pow(-1.,n+m-n+mp-n+1.);
          //n m -(mp+1)
          rv[n4+n2*m+mp+1-2*(mp+1-n)]           = rv[n4+n2*m+mp+1] *  pow(-1.,n+m-n+mp-n+1.);
          //n mp+1 -m
          rv[n4+n2*(mp+1)+m-2*(m-n)]            = rv[n4+n2*m+mp+1] *  pow(-1.,n+m-n+mp-n+1.);
          //n -(mp+1) m
          rv[n4+n2*(mp+1-2*(mp+1-n))+m]         = rv[n4+n2*m+mp+1] *  pow(-1.,n+m-n+mp-n+1.);
        }
      }
    }
  }
  //calculation of the coefficients for beta=acos(1/sqrt(3)) from the ones for beta=pi/2 via the flip transform. eq. A 24 in Gumerov 2012. Note, that there, different values for beta are claimed
  for(int n=0;n<truncation;n++) {
    for( int m=0;m<2*n+1;m++) {
      for(int mp=0;mp<2*n+1;mp++) {
        int n4 = (4*pow(n+0.,3.)-n)/3;
        int n2 = 2*n+1;
        rcoeff[n4+n2*m+mp]                      = 0;
        for(int y=0;y<(2*n)+1;y++) {
          rcoeff[n4+n2*m+mp]                   += rv[n4+n2*m+y]*rv[n4+n2*mp+y]* cos(acos(1/sqrt(3))*(y-n)+0.5*M_PI*(mp-n+m-n));
        }
      }
    }
  }
}


//calculates translation coefficients for a given level and q
void HSAXS::cal_transcoeff(int level, double qval, int tr) {
  //the structure of therecursion makes it necessary to calculate the translation coefficients to a higher truncation. See Gumerov, 2004, Chapters 3,7
  int p=2*tr-2;
  //length of the translation vector for a given target level (i.e. coarser level in the case of the upward pass
  double t=sqrt(0.75*pow(2.,-2.*level))*maxdist;
  tcoeff.resize((p+1)*(p+2)*(p+3)/6);
//analytical formula for m=n'=0. eq. A17 in Gumerov 2012
  for(int n=0;n<p+1;n++) {
    tcoeff[n*(n+1)*(n+2)/6]=pow(-1.,n+0.)*sqrt(2*n+1.)*gsl_sf_bessel_jl(n,qval*t);
  }
  //first recursion: advancement in m. eq. A18 in Gumerov 2012
  for(int m=0;m<tr;m++) {
    for(int n=(m+1);n<p;n++) {
      //the next three lines define offsets for n
      int i  = n*(n+1)*(n+2)/6;
      int ip = (n-1)*n*(n+1)/6;
      int in = (n+1)*(n+2)*(n+3)/6;
      //the next three line should define the offsets for m
      int j  = (m+1)*(n+1-0.5*m);
      int jp = m*(n-0.5*(m-1));
      int jn = m*(n+2-0.5*(m-1));
      tcoeff[i+j]  = 0;
      tcoeff[i+j] += 1/bv[(m+1)*(m+1)] * bv[n*n+m+n-2*m-1]     * tcoeff[ip+jp];
      tcoeff[i+j] -= 1/bv[(m+1)*(m+1)] * bv[(n+1)*(n+1)+m+n+1] * tcoeff[in+jn];
    }
  }
  //second recursion: advancement in np. eq. A19 in Gumerov 2012
  for(int m=0;m<tr;m++) {
    for(int np=m;np<p-1;np++) {
      for(int n=1;n<p;n++) {
        if(n>np) {
          //the next three lines define offsets for n
          int i  = n*(n+1)*(n+2)/6;
          int ip = (n-1)*n*(n+1)/6;
          int in = (n+1)*(n+2)*(n+3)/6;
          //the next three line should define the offsets for m
          int j  = m*(n+1-0.5*(m-1));
          int jp = m*(n-0.5*(m-1));
          int jn = m*(n+2-0.5*(m-1));
          tcoeff[i+j+np+1-m]              = 0;
          if(np>m) tcoeff[i+j+np+1-m]    += 1/av[np*np+m+np] * av[(np-1)*(np-1)+np+m-1] * tcoeff[i+j+np-1-m]; 
          tcoeff[i+j+np+1-m]             -= 1/av[np*np+m+np] * av[n*n+m+n]              * tcoeff[in+jn+np-(m)];
          if(n>m)  tcoeff[i+j+np+1-m]    += 1/av[np*np+m+np] * av[(n-1)*(n-1)+m+n-1]    * tcoeff[ip+jp+np-(m)];
        }
      }
    }
  }
}


//partial derivatives of the spherical basis functions
Vector2d HSAXS::dXHarmonics(unsigned k, unsigned int i, int n, int m, vector<Vector2d> &qRnm) {
  Vector2d                              dRdc =     (bv[n*(n+4)-m+1] * qRnm[p2*i+n*(n+2)+m+3] + bv[n*(n+2)+m+1] * qRnm[p2*i+n*(n+2)+m+1]);
  if((abs(m-n-1)<=(n-1))&&((n-1)>=0))   dRdc-=      bv[n*(n+2)-m]   * qRnm[p2*i+n*(n-2)+m-1];
  if((abs(m-n+1)<=(n-1))&&((n-1)>=0))   dRdc-=      bv[n*n+m]       * qRnm[p2*i+n*n-2*n+m+1];
  return dRdc;
}


Vector2d HSAXS::dYHarmonics(unsigned k, unsigned int i, int n, int m, vector<Vector2d> &qRnm) {
  Vector2d                              dRdc =     (bv[n*(n+4)-m+1] * qRnm[p2*i+n*(n+2)+m+3] - bv[n*(n+2)+m+1] * qRnm[p2*i+n*(n+2)+m+1]);
  if((abs(m-n-1)<=(n-1))&&((n-1)>=0))   dRdc+=      bv[n*(n+2)-m]   * qRnm[p2*i+n*(n-2)+m-1];
  if((abs(m-n+1)<=(n-1))&&((n-1)>=0))   dRdc-=      bv[n*n+m]       * qRnm[p2*i+n*(n-2)+m+1];
  return dRdc;
}


Vector2d HSAXS::dZHarmonics(unsigned k, unsigned int i, int n, int m, vector<Vector2d> &qRnm) {
  Vector2d                              dRdc =      av[n*n+m]       * qRnm[p2*i+n*(n+2)+m+2];
  if((abs(m-n)<=(n-1))&&((n-1)>=0))     dRdc-=      av[n*(n-2)+m]   * qRnm[p2*i+n*(n-2)+m];
  return dRdc;
}


//below, a bunch a functions for finding the indices/centres of a box or parent box follow. Mostly based on bit(de)interleaving and bitshifts

//calculates a binary representation for floating pont decimals of format 0.x. 
std::bitset<24> HSAXS::cal_bin(double dec) {
  std::bitset<24> Bin;
  for(unsigned j=0;j<24;j++) {
    dec *=2;
    if(dec>=1.) { 
      dec-=1;
      Bin[23-j]=1;
    }
    else {
      Bin[23-j]=0;
    }
  }
  return Bin;
}


//use bit interleaving to combine three coordinate indices to one index
std::bitset<72> HSAXS::interleave(std::bitset<24> x, std::bitset<24> y, std::bitset<24> z) {
  std::bitset<72> xyz;
  for(unsigned j=0;j<24;j++) {
    xyz[3*(23-j)+2] = x[23-j];
    xyz[3*(23-j)+1] = y[23-j];
    xyz[3*(23-j)+0] = z[23-j];
  }
  return xyz;
}


//revert bit interleaving
void HSAXS::deinterleave(std::bitset<72> xyz, std::bitset<24> &x, std::bitset<24> &y, std::bitset<24> &z) {
  for(unsigned j=0;j<24;j++) {
    x[23-j] = xyz[3*(23-j)+2];
    y[23-j] = xyz[3*(23-j)+1];
    z[23-j] = xyz[3*(23-j)+0];
  }
}


//calculate the box centre for a given box index at level l
Vector HSAXS::cal_boxcentre(unsigned index, unsigned l) {
  std::bitset<24> n(index);
  Vector coord;
  std::bitset<8> nx;
  std::bitset<8> ny;
  std::bitset<8> nz;
  for(unsigned j=0;j<8;j++) {
    nx[7-j] = n[3*(7-j)+2];
    ny[7-j] = n[3*(7-j)+1];
    nz[7-j] = n[3*(7-j)+0];
  }
  nx <<=(1);
  ny <<=(1);
  nz <<=(1);
  nx[0]=1;
  ny[0]=1;
  nz[0]=1;
  for(int j=l;j>=0;j--) {
    coord[0] += nx[j]*pow(2.,-(l-j+1.));
    coord[1] += ny[j]*pow(2.,-(l-j+1.));
    coord[2] += nz[j]*pow(2.,-(l-j+1.));
  }
  return coord;
}


//calculate the integer box index at level l pertaining to an atom with interleaved binary coordinates coord
long unsigned HSAXS::l_index(unsigned l, std::bitset<72> coord) {
  coord >>=(72-3*l);
  return coord.to_ulong();//index;
}


//calculate integer box index for an atom with coordinates x,y,z at level l
long unsigned HSAXS::own_index(unsigned ownLevel, double x, double y, double z) {
  return l_index(ownLevel,interleave(cal_bin(x),cal_bin(y),cal_bin(z)));
}


//calculate integer box index for the parent box of an atom with coordinates x,y,z at level l
long unsigned HSAXS::parent_index(long unsigned index) {
  std::bitset<24> n(index);
  n >>=3;
  return n.to_ulong();
}


long unsigned HSAXS::parent_index(unsigned ownLevel, double x, double y, double z) {
  long unsigned index=own_index(ownLevel, x, y, z);
  std::bitset<24> n(index);
  n >>=3;
  return n.to_ulong();
}


//calculate box centre for the parent box of an atom with coordinates x,y,z at level l
Vector HSAXS::parent_centre(unsigned ownLevel, double x, double y, double z) {
  return cal_boxcentre(parent_index(own_index(ownLevel, x, y, z)), ownLevel-1);
}


//calculate box centre for an atom with coordinates x,y,z at level l
Vector HSAXS::own_centre(unsigned ownLevel, double x, double y, double z) {
  return cal_boxcentre(own_index(ownLevel, x, y, z), ownLevel);
}


//formfactors for MARTINI calculations
void HSAXS::getMartiniSFparam(const vector<AtomNumber> &atoms, vector<vector<long double> > &parameter)
{
  vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
  if( moldat.size()==1 ){
    log<<"  MOLINFO DATA found, using proper atom names\n";
    for(unsigned i=0;i<atoms.size();++i){
      string Aname = moldat[0]->getAtomName(atoms[i]);
      string Rname = moldat[0]->getResidueName(atoms[i]);
      if(Rname=="ALA") {
        if(Aname=="BB") {
          parameter[i].push_back(9.045);
          parameter[i].push_back(-0.098114);
          parameter[i].push_back(7.54281);
          parameter[i].push_back(-1.97438);
          parameter[i].push_back(-8.32689);
          parameter[i].push_back(6.09318);
          parameter[i].push_back(-1.18913);
        } else error("Atom name not known");
      } else if(Rname=="ARG") {
        if(Aname=="BB") {
          parameter[i].push_back(10.729);
          parameter[i].push_back(-0.0392574);
          parameter[i].push_back(1.15382);
          parameter[i].push_back(-0.155999);
          parameter[i].push_back(-2.43619);
          parameter[i].push_back(1.72922);
          parameter[i].push_back(-0.33799);
        } else if(Aname=="SC1"){
          parameter[i].push_back(-2.797);
          parameter[i].push_back(0.472403);
          parameter[i].push_back(8.07424);
          parameter[i].push_back(4.37299);
          parameter[i].push_back(-10.7398);
          parameter[i].push_back(4.95677);
          parameter[i].push_back(-0.725797);
        } else if(Aname=="SC2"){
          parameter[i].push_back(15.396);
          parameter[i].push_back(0.0636736);
          parameter[i].push_back(-1.258);
          parameter[i].push_back(1.93135);
          parameter[i].push_back(-4.45031);
          parameter[i].push_back(2.49356);
          parameter[i].push_back(-0.410721);
        } else error("Atom name not known");
      } else if(Rname=="ASN") {
        if(Aname=="BB") {
          parameter[i].push_back(10.738);
          parameter[i].push_back(-0.0402162);
          parameter[i].push_back(1.03007);
          parameter[i].push_back(-0.254174);
          parameter[i].push_back(-2.12015);
          parameter[i].push_back(1.55535);
          parameter[i].push_back(-0.30963);
        } else if(Aname=="SC1"){
          parameter[i].push_back(9.249);
          parameter[i].push_back(-0.0148678);
          parameter[i].push_back(5.52169);
          parameter[i].push_back(0.00853212);
          parameter[i].push_back(-6.71992);
          parameter[i].push_back(3.93622);
          parameter[i].push_back(-0.64973);
        } else error("Atom name not known");
      } else if(Rname=="ASP") {
        if(Aname=="BB") {
          parameter[i].push_back(10.695);
          parameter[i].push_back(-0.0410247);
          parameter[i].push_back(1.03656);
          parameter[i].push_back(-0.298558);
          parameter[i].push_back(-2.06064);
          parameter[i].push_back(1.53495);
          parameter[i].push_back(-0.308365);
        } else if(Aname=="SC1"){
          parameter[i].push_back(9.476);
          parameter[i].push_back(-0.0254664);
          parameter[i].push_back(5.57899);
          parameter[i].push_back(-0.395027);
          parameter[i].push_back(-5.9407);
          parameter[i].push_back(3.48836);
          parameter[i].push_back(-0.569402);
        } else error("Atom name not known");
      } else if(Rname=="CYS") {
        if(Aname=="BB") {
          parameter[i].push_back(10.698);
          parameter[i].push_back(-0.0233493);
          parameter[i].push_back(1.18257);
          parameter[i].push_back(0.0684463);
          parameter[i].push_back(-2.792);
          parameter[i].push_back(1.88995);
          parameter[i].push_back(-0.360229);
        } else if(Aname=="SC1"){
          parameter[i].push_back(8.199);
          parameter[i].push_back(-0.0261569);
          parameter[i].push_back(6.79677);
          parameter[i].push_back(-0.343845);
          parameter[i].push_back(-5.03578);
          parameter[i].push_back(2.7076);
          parameter[i].push_back(-0.420714);
        } else error("Atom name not known");
      } else if(Rname=="GLN") {
        if(Aname=="BB") {
          parameter[i].push_back(10.728);
          parameter[i].push_back(-0.0391984);
          parameter[i].push_back(1.09264);
          parameter[i].push_back(-0.261555);
          parameter[i].push_back(-2.21245);
          parameter[i].push_back(1.62071);
          parameter[i].push_back(-0.322325);
        } else if(Aname=="SC1"){
          parameter[i].push_back(8.317);
          parameter[i].push_back(-0.229045);
          parameter[i].push_back(12.6338);
          parameter[i].push_back(-7.6719);
          parameter[i].push_back(-5.8376);
          parameter[i].push_back(5.53784);
          parameter[i].push_back(-1.12604);
        } else error("Atom name not known");
      } else if(Rname=="GLU") {
        if(Aname=="BB") {
          parameter[i].push_back(10.694);
          parameter[i].push_back(-0.0521961);
          parameter[i].push_back(1.11153);
          parameter[i].push_back(-0.491995);
          parameter[i].push_back(-1.86236);
          parameter[i].push_back(1.45332);
          parameter[i].push_back(-0.29708);
        } else if(Aname=="SC1"){
          parameter[i].push_back(8.544);
          parameter[i].push_back(-0.249555);
          parameter[i].push_back(12.8031);
          parameter[i].push_back(-8.42696);
          parameter[i].push_back(-4.66486);
          parameter[i].push_back(4.90004);
          parameter[i].push_back(-1.01204);
        } else error("Atom name not known");
      } else if(Rname=="GLY") {
        if(Aname=="BB") {
          parameter[i].push_back(9.977);
          parameter[i].push_back(-0.0285799);
          parameter[i].push_back(1.84236);
          parameter[i].push_back(-0.0315192);
          parameter[i].push_back(-2.88326);
          parameter[i].push_back(1.87323);
          parameter[i].push_back(-0.345773);
        } else error("Atom name not known");
      } else if(Rname=="HIS") {
        if(Aname=="BB") {
          parameter[i].push_back(10.721);
          parameter[i].push_back(-0.0379337);
          parameter[i].push_back(1.06028);
          parameter[i].push_back(-0.236143);
          parameter[i].push_back(-2.17819);
          parameter[i].push_back(1.58357);
          parameter[i].push_back(-0.31345);
        } else if(Aname=="SC1"){
          parameter[i].push_back(-0.424);
          parameter[i].push_back(0.665176);
          parameter[i].push_back(3.4369);
          parameter[i].push_back(2.93795);
          parameter[i].push_back(-5.18288);
          parameter[i].push_back(2.12381);
          parameter[i].push_back(-0.284224);
        } else if(Aname=="SC2"){
          parameter[i].push_back(5.363);
          parameter[i].push_back(-0.0176945);
          parameter[i].push_back(2.9506);
          parameter[i].push_back(-0.387018);
          parameter[i].push_back(-1.83951);
          parameter[i].push_back(0.9703);
          parameter[i].push_back(-0.1458);
        } else if(Aname=="SC3"){
          parameter[i].push_back(5.784);
          parameter[i].push_back(-0.0293129);
          parameter[i].push_back(2.74167);
          parameter[i].push_back(-0.520875);
          parameter[i].push_back(-1.62949);
          parameter[i].push_back(0.902379);
          parameter[i].push_back(-0.139957);
        } else error("Atom name not known");
      } else if(Rname=="ILE") {
        if(Aname=="BB") {
          parameter[i].push_back(10.699);
          parameter[i].push_back(-0.0188962);
          parameter[i].push_back(1.217);
          parameter[i].push_back(0.242481);
          parameter[i].push_back(-3.13898);
          parameter[i].push_back(2.07916);
          parameter[i].push_back(-0.392574);
        } else if(Aname=="SC1"){
          parameter[i].push_back(-4.448);
          parameter[i].push_back(1.20996);
          parameter[i].push_back(11.5141);
          parameter[i].push_back(6.98895);
          parameter[i].push_back(-19.1948);
          parameter[i].push_back(9.89207);
          parameter[i].push_back(-1.60877);
        } else error("Atom name not known");
      } else if(Rname=="LEU") {
        if(Aname=="BB") {
          parameter[i].push_back(10.692);
          parameter[i].push_back(-0.209448);
          parameter[i].push_back(1.73738);
          parameter[i].push_back(-1.33726);
          parameter[i].push_back(-1.3065);
          parameter[i].push_back(1.25273);
          parameter[i].push_back(-0.265001);
        } else if(Aname=="SC1"){
          parameter[i].push_back(-4.448);
          parameter[i].push_back(2.1063);
          parameter[i].push_back(6.72381);
          parameter[i].push_back(14.6954);
          parameter[i].push_back(-23.7197);
          parameter[i].push_back(10.7247);
          parameter[i].push_back(-1.59146);
        } else error("Atom name not known");
      } else if(Rname=="LYS") {
        if(Aname=="BB") {
          parameter[i].push_back(10.706);
          parameter[i].push_back(-0.0468629);
          parameter[i].push_back(1.09477);
          parameter[i].push_back(-0.432751);
          parameter[i].push_back(-1.94335);
          parameter[i].push_back(1.49109);
          parameter[i].push_back(-0.302589);
        } else if(Aname=="SC1"){
          parameter[i].push_back(-2.796);
          parameter[i].push_back(0.508044);
          parameter[i].push_back(7.91436);
          parameter[i].push_back(4.54097);
          parameter[i].push_back(-10.8051);
          parameter[i].push_back(4.96204);
          parameter[i].push_back(-0.724414);
        } else if(Aname=="SC2"){
          parameter[i].push_back(3.070);
          parameter[i].push_back(-0.0101448);
          parameter[i].push_back(4.67994);
          parameter[i].push_back(-0.792529);
          parameter[i].push_back(-2.09142);
          parameter[i].push_back(1.02933);
          parameter[i].push_back(-0.137787);
        } else error("Atom name not known");
      } else if(Rname=="MET") {
        if(Aname=="BB") {
          parameter[i].push_back(10.671);
          parameter[i].push_back(-0.0433724);
          parameter[i].push_back(1.13784);
          parameter[i].push_back(-0.40768);
          parameter[i].push_back(-2.00555);
          parameter[i].push_back(1.51673);
          parameter[i].push_back(-0.305547);
        } else if(Aname=="SC1"){
          parameter[i].push_back(5.85);
          parameter[i].push_back(-0.0485798);
          parameter[i].push_back(17.0391);
          parameter[i].push_back(-3.65327);
          parameter[i].push_back(-13.174);
          parameter[i].push_back(8.68286);
          parameter[i].push_back(-1.56095);
        } else error("Atom name not known");
      } else if(Rname=="PHE") {
        if(Aname=="BB") {
          parameter[i].push_back(10.741);
          parameter[i].push_back(-0.0317276);
          parameter[i].push_back(1.15599);
          parameter[i].push_back(0.0276186);
          parameter[i].push_back(-2.74757);
          parameter[i].push_back(1.88783);
          parameter[i].push_back(-0.363525);
        } else if(Aname=="SC1"){
          parameter[i].push_back(-0.636);
          parameter[i].push_back(0.527882);
          parameter[i].push_back(6.77612);
          parameter[i].push_back(3.18508);
          parameter[i].push_back(-8.92826);
          parameter[i].push_back(4.29752);
          parameter[i].push_back(-0.65187);
        } else if(Aname=="SC2"){
          parameter[i].push_back(-0.424);
          parameter[i].push_back(0.389174);
          parameter[i].push_back(4.11761);
          parameter[i].push_back(2.29527);
          parameter[i].push_back(-4.7652);
          parameter[i].push_back(1.97023);
          parameter[i].push_back(-0.262318);
        } else if(Aname=="SC3"){
          parameter[i].push_back(-0.424);
          parameter[i].push_back(0.38927);
          parameter[i].push_back(4.11708);
          parameter[i].push_back(2.29623);
          parameter[i].push_back(-4.76592);
          parameter[i].push_back(1.97055);
          parameter[i].push_back(-0.26238);
        } else error("Atom name not known");
      } else if(Rname=="PRO") {
        if(Aname=="BB") {
          parameter[i].push_back(11.434);
          parameter[i].push_back(-0.033323);
          parameter[i].push_back(0.472014);
          parameter[i].push_back(-0.290854);
          parameter[i].push_back(-1.81409);
          parameter[i].push_back(1.39751);
          parameter[i].push_back(-0.280407);
        } else if(Aname=="SC1"){
          parameter[i].push_back(-2.796);
          parameter[i].push_back(0.95668);
          parameter[i].push_back(6.84197);
          parameter[i].push_back(6.43774);
          parameter[i].push_back(-12.5068);
          parameter[i].push_back(5.64597);
          parameter[i].push_back(-0.825206);
        } else error("Atom name not known");
      } else if(Rname=="SER") {
        if(Aname=="BB") {
          parameter[i].push_back(10.699);
          parameter[i].push_back(-0.0325828);
          parameter[i].push_back(1.20329);
          parameter[i].push_back(-0.0674351);
          parameter[i].push_back(-2.60749);
          parameter[i].push_back(1.80318);
          parameter[i].push_back(-0.346803);
        } else if(Aname=="SC1"){
          parameter[i].push_back(3.298);
          parameter[i].push_back(-0.0366801);
          parameter[i].push_back(5.11077);
          parameter[i].push_back(-1.46774);
          parameter[i].push_back(-1.48421);
          parameter[i].push_back(0.800326);
          parameter[i].push_back(-0.108314);
        } else error("Atom name not known");
      } else if(Rname=="THR") {
        if(Aname=="BB") {
          parameter[i].push_back(10.697);
          parameter[i].push_back(-0.0242955);
          parameter[i].push_back(1.24671);
          parameter[i].push_back(0.146423);
          parameter[i].push_back(-2.97429);
          parameter[i].push_back(1.97513);
          parameter[i].push_back(-0.371479);
        } else if(Aname=="SC1"){
          parameter[i].push_back(2.366);
          parameter[i].push_back(0.0297604);
          parameter[i].push_back(11.9216);
          parameter[i].push_back(-9.32503);
          parameter[i].push_back(1.9396);
          parameter[i].push_back(0.0804861);
          parameter[i].push_back(-0.0302721);
        } else error("Atom name not known");
      } else if(Rname=="TRP") {
        if(Aname=="BB") {
          parameter[i].push_back(10.689);
          parameter[i].push_back(-0.0265879);
          parameter[i].push_back(1.17819);
          parameter[i].push_back(0.0386457);
          parameter[i].push_back(-2.75634);
          parameter[i].push_back(1.88065);
          parameter[i].push_back(-0.360217);
        } else if(Aname=="SC1"){
          parameter[i].push_back(0.084);
          parameter[i].push_back(0.752407);
          parameter[i].push_back(5.3802);
          parameter[i].push_back(4.09281);
          parameter[i].push_back(-9.28029);
          parameter[i].push_back(4.45923);
          parameter[i].push_back(-0.689008);
        } else if(Aname=="SC2"){
          parameter[i].push_back(5.739);
          parameter[i].push_back(0.0298492);
          parameter[i].push_back(4.60446);
          parameter[i].push_back(1.34463);
          parameter[i].push_back(-5.69968);
          parameter[i].push_back(2.84924);
          parameter[i].push_back(-0.433781);
        } else if(Aname=="SC3"){
          parameter[i].push_back(-0.424);
          parameter[i].push_back(0.388576);
          parameter[i].push_back(4.11859);
          parameter[i].push_back(2.29485);
          parameter[i].push_back(-4.76255);
          parameter[i].push_back(1.96849);
          parameter[i].push_back(-0.262015);
        } else if(Aname=="SC4"){
          parameter[i].push_back(-0.424);
          parameter[i].push_back(0.387685);
          parameter[i].push_back(4.12153);
          parameter[i].push_back(2.29144);
          parameter[i].push_back(-4.7589);
          parameter[i].push_back(1.96686);
          parameter[i].push_back(-0.261786);
        } else error("Atom name not known");
      } else if(Rname=="TYR") {
        if(Aname=="BB") {
          parameter[i].push_back(10.689);
          parameter[i].push_back(-0.0193526);
          parameter[i].push_back(1.18241);
          parameter[i].push_back(0.207318);
          parameter[i].push_back(-3.0041);
          parameter[i].push_back(1.99335);
          parameter[i].push_back(-0.376482);
        } else if(Aname=="SC1"){
          parameter[i].push_back(-0.636);
          parameter[i].push_back(0.528902);
          parameter[i].push_back(6.78168);
          parameter[i].push_back(3.17769);
          parameter[i].push_back(-8.93667);
          parameter[i].push_back(4.30692);
          parameter[i].push_back(-0.653993);
        } else if(Aname=="SC2"){
          parameter[i].push_back(-0.424);
          parameter[i].push_back(0.388811);
          parameter[i].push_back(4.11851);
          parameter[i].push_back(2.29545);
          parameter[i].push_back(-4.7668);
          parameter[i].push_back(1.97131);
          parameter[i].push_back(-0.262534);
        } else if(Aname=="SC3"){
          parameter[i].push_back(4.526);
          parameter[i].push_back(-0.00381305);
          parameter[i].push_back(5.8567);
          parameter[i].push_back(-0.214086);
          parameter[i].push_back(-4.63649);
          parameter[i].push_back(2.52869);
          parameter[i].push_back(-0.39894);
        } else error("Atom name not known");
      } else if(Rname=="VAL") {
        if(Aname=="BB") {
          parameter[i].push_back(10.691);
          parameter[i].push_back(-0.0162929);
          parameter[i].push_back(1.24446);
          parameter[i].push_back(0.307914);
          parameter[i].push_back(-3.27446);
          parameter[i].push_back(2.14788);
          parameter[i].push_back(-0.403259);
        } else if(Aname=="SC1"){
          parameter[i].push_back(-3.516);
          parameter[i].push_back(1.62307);
          parameter[i].push_back(5.43064);
          parameter[i].push_back(9.28809);
          parameter[i].push_back(-14.9927);
          parameter[i].push_back(6.6133);
          parameter[i].push_back(-0.964977);
        } else error("Atom name not known");
      } else error("Residue not known");
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
}
//atomistic formfactors
void HSAXS::calculateASF(const vector<AtomNumber> &atoms, vector<vector<long double> > &FF_tmp)
{
   enum { H, C, N, O, P, S, X, NTT };
   map<string, unsigned> AA_map;
   AA_map["H"] = H;
   AA_map["C"] = C;
   AA_map["N"] = N;
   AA_map["O"] = O;
   AA_map["P"] = P;
   AA_map["S"] = S;
   AA_map["X"] = X;
 
   vector<vector<double> > param_a;
   vector<vector<double> > param_b;
   vector<double> param_c;
   vector<double> param_v;

   param_a.resize(NTT, vector<double>(5));
   param_b.resize(NTT, vector<double>(5));
   param_c.resize(NTT);
   param_v.resize(NTT);

   param_a[H][0] = 0.493002; param_b[H][0] = 10.5109; param_c[H] = 0.003038;
   param_a[H][1] = 0.322912; param_b[H][1] = 26.1257; param_v[H] = 5.15;
   param_a[H][2] = 0.140191; param_b[H][2] = 3.14236;
   param_a[H][3] = 0.040810; param_b[H][3] = 57.7997;
   param_a[H][4] = 0.0;      param_b[H][4] = 1.0;

   param_a[C][0] = 2.31000; param_b[C][0] = 20.8439; param_c[C] = 0.215600;
   param_a[C][1] = 1.02000; param_b[C][1] = 10.2075; param_v[C] = 16.44;
   param_a[C][2] = 1.58860; param_b[C][2] = 0.56870;
   param_a[C][3] = 0.86500; param_b[C][3] = 51.6512;
   param_a[C][4] = 0.0;     param_b[C][4] = 1.0;

   param_a[N][0] = 12.2126; param_b[N][0] = 0.00570; param_c[N] = -11.529;
   param_a[N][1] = 3.13220; param_b[N][1] = 9.89330; param_v[N] = 2.49;
   param_a[N][2] = 2.01250; param_b[N][2] = 28.9975;
   param_a[N][3] = 1.16630; param_b[N][3] = 0.58260;
   param_a[N][4] = 0.0;     param_b[N][4] = 1.0;

   param_a[O][0] = 3.04850; param_b[O][0] = 13.2771; param_c[O] = 0.250800 ;
   param_a[O][1] = 2.28680; param_b[O][1] = 5.70110; param_v[O] = 9.13;
   param_a[O][2] = 1.54630; param_b[O][2] = 0.32390;
   param_a[O][3] = 0.86700; param_b[O][3] = 32.9089;
   param_a[O][4] = 0.0;     param_b[O][4] = 1.0;

   param_a[P][0] = 6.43450; param_b[P][0] = 1.90670; param_c[P] = 1.11490;
   param_a[P][1] = 4.17910; param_b[P][1] = 27.1570; param_v[P] = 5.73;
   param_a[P][2] = 1.78000; param_b[P][2] = 0.52600;
   param_a[P][3] = 1.49080; param_b[P][3] = 68.1645;
   param_a[P][4] = 0.0;     param_b[P][4] = 1.0;

   param_a[S][0] = 6.90530; param_b[S][0] = 1.46790; param_c[S] = 0.866900;
   param_a[S][1] = 5.20340; param_b[S][1] = 22.2151; param_v[S] = 19.86;
   param_a[S][2] = 1.43790; param_b[S][2] = 0.25360;
   param_a[S][3] = 1.58630; param_b[S][3] = 56.1720;
   param_a[S][4] = 0.0;     param_b[S][4] = 1.0;

   param_a[X][0] = 0; param_b[X][0] = 0; param_c[X] = 0;
   param_a[X][1] = 0; param_b[X][1] = 0; param_v[X] = 0;
   param_a[X][2] = 0; param_b[X][2] = 0;
   param_a[X][3] = 0; param_b[X][3] = 0;
   param_a[X][4] = 0; param_b[X][4] = 0;

   vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();

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
       std::string type_s = std::string(1,type);
       if(AA_map.find(type_s) != AA_map.end()){
         const unsigned index=AA_map[type_s];
         const double rho = 0.334;
         const double volr = pow(param_v[index], (2.0/3.0)) /(4. * M_PI);
         for(unsigned k=0;k<numq;++k){
           const double q = q_list[k];
           const double s = q / (4. * M_PI);
           FF_tmp[k][i] = param_c[index];
           // SUM [a_i * EXP( - b_i * (q/4pi)^2 )] Waasmaier and Kirfel (1995)
           for(unsigned j=0;j<5;j++) {
             FF_tmp[k][i] += param_a[index][j]*exp(-param_b[index][j]*s*s);
           }
           // subtract solvation: rho * v_i * EXP( (- v_i^(2/3) / (4pi)) * q^2  )
           FF_tmp[k][i] -= rho*param_v[index]*exp(-volr*q*q);
         }
       } else {
         error("Wrong atom type "+type_s+" from atom name "+name+"\n"); 
       }
     }
   } else {
     error("MOLINFO DATA not found\n");
   }
}

}
}
