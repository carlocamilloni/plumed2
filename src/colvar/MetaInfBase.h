/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#ifndef __PLUMED_colvar_MetaInfBase_h
#define __PLUMED_colvar_MetaInfBase_h

#include "core/Atoms.h"
#include "tools/Random.h"
#include "tools/File.h"

#include <vector>
#include <string>

using namespace std;

namespace PLMD {
    class Keywords;

    class MetaInfBase {
        private:
            const double sqrt2_div_pi;

            // experimental values
            vector<double> parameters;

            // noise type
            unsigned noise_type_;
            enum { GAUSS, MGAUSS, OUTLIERS };

            // scale is data scaling factor
            unsigned scale_prior_;
            enum { SC_GAUSS, SC_FLAT };
            bool   doscale_;
            double scale_;
            double scale_mu_;
            double scale_sigma_;
            double scale_min_;
            double scale_max_;
            double Dscale_;

            // sigma is data uncertainty
            vector<double> sigma_;
            double sigma_min_;
            double sigma_max_;
            double Dsigma_;

            // sigma_mean is uncertainty in the mean estimate
            vector<double> sigma_mean_;
            vector<double> variance_;

            // sigma_mean rescue params
            double sm_mod_;
            double sm_mod_min_;
            double sm_mod_max_;
            double Dsm_mod_;
            double max_force_;
            double fmax_;

            // temperature in kbt
            double   kbt_;

            // number of data points
            unsigned ndata_;

            // Monte Carlo stuff
            vector<Random> random;
            unsigned MCsteps_;
            unsigned MCstride_;
            long unsigned MCaccept_;
            long unsigned MCtrial_;
            double accept;

            // restart
            unsigned write_stride_;

            // others
            bool     master;
            bool     do_reweight;
            bool     do_optsigmamean_;
            bool     do_mc_single_;
            unsigned nrep_;
            unsigned replica_;
            unsigned narg;
            vector<double> output_force;

            double getEnergySPE(const vector<double> &mean,
                                const vector<double> &sigma,
                                const double scale);
            double getEnergyGJE(const vector<double> &mean,
                                const vector<double> &sigma,
                                const double scale);
            void doMonteCarlo(const vector<double> &mean_,
                              Communicator& comm,
                              Communicator& multi_sim_comm);
            double getEnergyForceSPE(const vector<double> &arguments,
                                     Communicator& comm,
                                     Communicator& multi_sim_comm,
                                     const vector<double> &mean,
                                     const double fact);
            double getEnergyForceGJE(const vector<double> &arguments,
                                     Communicator& comm,
                                     Communicator& multi_sim_comm,
                                     const vector<double> &mean,
                                     const double fact);
            void writeStatus(double timestep, unsigned step, OFile& sfile_);
            void optSigmaMean(Communicator& comm,
                              Communicator& multi_sim_comm);
        
        public:
            static void registerKeywords(Keywords& keys);
            MetaInfBase();
            ~MetaInfBase();
            vector<double>& getOutputForce();
            void set(const std::string& definition,
                     std::string& errormsg,
                     vector<double>& datapoints,
                     bool restart,
                     double kBoltzmann,
                     double kbt,
                     Communicator& comm,
                     Communicator& multi_sim_comm,
                     Log& log,
                     const std::string& label,
                     IFile& restart_sfile,
                     OFile& sfile);
            double calculate(vector<double>& arguments,
                             double timestep,
                             const long int step,
                             const bool exchange_step,
                             const bool checkpoint,
                             Communicator& comm,
                             Communicator& multi_sim_comm,
                             OFile& sfile);
    };
}

#endif
