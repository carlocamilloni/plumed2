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

/* TODO: multiple timestep [x] */
/*       statusfiles       [x] */
/*       sigma_mean opt    [ ] */
/*       bibliography      [ ] */
/*       random numbers    [x] */
/*       other components  [ ] */

#include "MetaInfBase.h"
#include "tools/Keywords.h"

#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "core/Value.h"
#include "tools/File.h"
#include "tools/Random.h"
#include <cmath>
#include <ctime>

using namespace std;

namespace PLMD {
    void MetaInfBase::registerKeywords(Keywords& keys) {
        keys.addFlag("REWEIGHT",false,"simple REWEIGHT using the latest ARG as energy"); 
        keys.addFlag("SCALEDATA",false,"Set to TRUE if you want to sample a scaling factor common to all values and replicas");  
        keys.addFlag("OPTSIGMAMEAN",false,"Set to TRUE if you want to scale sigma_mean to prevent too high forces");  
        keys.addFlag("MCSINGLE",false,"Set to TRUE if you want to change a single sigma per MC move (only for NOISETYPE=MGAUSS)");  
        keys.add("compulsory","NOISETYPE","functional form of the noise (GAUSS,MGAUSS,OUTLIERS)");
        keys.add("compulsory","SCALE_PRIOR","FLAT","either FLAT or GAUSSIAN");
        keys.add("compulsory","SCALE0","initial value of the uncertainty parameter");
        keys.add("optional","SCALE_MIN","minimum value of the uncertainty parameter");
        keys.add("optional","SCALE_MAX","maximum value of the uncertainty parameter");
        keys.add("optional","SCALE_SIGMA","maximum value of the uncertainty parameter");
        keys.add("optional","DSCALE","maximum MC move of the uncertainty parameter");
        keys.add("compulsory","SIGMA0","initial value of the uncertainty parameter");
        keys.add("compulsory","SIGMA_MIN","minimum value of the uncertainty parameter");
        keys.add("compulsory","SIGMA_MAX","maximum value of the uncertainty parameter");
        keys.add("optional","DSIGMA","maximum MC move of the uncertainty parameter");
        keys.add("optional","SIGMA_MEAN0","starting value for the uncertainty in the mean estimate");
        keys.add("optional","SIGMA_MEAN_MOD0","starting value for sm modifier");
        keys.add("optional","SIGMA_MEAN_MOD_MIN","starting value for sm modifier");
        keys.add("optional","SIGMA_MEAN_MOD_MAX","starting value for sm modifier");
        keys.add("optional","DSIGMA_MEAN_MOD","step value for sm modifier");
        keys.add("optional","MAX_FORCE","maximum allowable force");
        keys.add("optional","TEMP","the system temperature - this is only needed if code doesnt' pass the temperature to plumed");
        keys.add("optional","MC_STEPS","number of MC steps");
        keys.add("optional","MC_STRIDE","MC stride");
        keys.add("optional","STRIDE","calculation stride");
        keys.add("optional","STATUS_FILE","write a file with all the data usefull for restart/continuation of Metainference");
        keys.add("compulsory","WRITE_STRIDE","write the status to a file every N steps, this can be used for restart/continuation");
    }

    void MetaInfBase::set(const std::string& definition,
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
                          OFile& sfile) {

        parameters = datapoints;
        narg = datapoints.size();
        ndata_ = datapoints.size();
        output_force.resize(ndata_);

        vector<string> data = Tools::getWords(definition);
        if (data.size() < 1) {
            errormsg = "Missing all input for MetaInfBase";
        }

        // set up replica stuff 
        master = (comm.Get_rank() == 0);
        if (master) {
            nrep_    = multi_sim_comm.Get_size();
            replica_ = multi_sim_comm.Get_rank();
        } else {
            nrep_    = 0;
            replica_ = 0;
        }
        comm.Sum(&nrep_, 1);
        comm.Sum(&replica_, 1);

        // Start the actual parsing
        Tools::parseFlag(data, "REWEIGHT", do_reweight);

        string stringa_noise;
        Tools::parse(data, "NOISETYPE", stringa_noise);
        if (stringa_noise=="GAUSS") {
            noise_type_ = GAUSS;
        } else if (stringa_noise == "MGAUSS") {
            noise_type_ = MGAUSS;
        } else if (stringa_noise == "OUTLIERS") {
            noise_type_ = OUTLIERS;
        } else {
            errormsg = "Unknown noise type!"; 
        }

        // Status file
        Tools::parse(data, "WRITE_STRIDE", write_stride_);
        string status_file_name_;
        Tools::parse(data, "STATUS_FILE", status_file_name_);
        if (status_file_name_ == "") {
            status_file_name_ = "MISTATUS" + label;
        }

        Tools::parseFlag(data, "OPTSIGMAMEAN", do_optsigmamean_);
        Tools::parseFlag(data, "MCSINGLE", do_mc_single_);
        if (noise_type_ != MGAUSS && do_mc_single_) {
            errormsg = "MCSINGLE is only available with MGAUSS";
        }

        Tools::parseFlag(data, "SCALEDATA", doscale_);
        if (doscale_) {
            string stringa_noise;
            Tools::parse(data, "SCALE_PRIOR", stringa_noise);
            if (stringa_noise == "GAUSSIAN") {
                scale_prior_ = SC_GAUSS; 
            } else if (stringa_noise == "FLAT") {
                scale_prior_ = SC_FLAT;
            } else {
                errormsg = "Unknown SCALE_PRIOR type!";
            }

            Tools::parse(data, "SCALE0", scale_);
            if(scale_prior_ == SC_GAUSS) {
                Tools::parse(data, "SCALE_SIGMA", scale_sigma_);
                scale_mu_ = scale_;
            } else {
                Tools::parse(data, "SCALE_MIN", scale_min_);
                Tools::parse(data, "SCALE_MAX", scale_max_);
            }

            if (scale_prior_ == SC_GAUSS && scale_sigma_ < 0.) {
                errormsg = "SCALE_SIGMA must be set when using SCALE_PRIOR=GAUSS";
            }
            if (scale_prior_ == SC_FLAT && scale_max_ < scale_min_) {
                errormsg = "SCALE_MAX and SCALE_MIN must be set when using SCALE_PRIOR=FLAT";
            }

            Tools::parse(data, "DSCALE", Dscale_);
            if (Dscale_ < 0) {
                if (scale_prior_ == SC_FLAT) {
                    Dscale_ = 0.05 * (scale_max_ - scale_min_);
                }
                if (scale_prior_ == SC_GAUSS) {
                    Dscale_ = 0.05 * scale_sigma_;
                }
            }
        } else {
            scale_ = 1.0;
        }

        vector<double> readsigma;
        Tools::parseVector(data, "SIGMA0", readsigma);
        if (noise_type_ != MGAUSS && readsigma.size() > 1) {
            errormsg = "If you want to use more than one SIGMA you should use NOISETYPE=MGAUSS";
        }
        if (noise_type_ == MGAUSS) {
            if (readsigma.size() == narg) {
                sigma_.resize(narg);
                sigma_ = readsigma;
            } else if (readsigma.size() == 1) {
                sigma_.resize(narg, readsigma[0]);
            } else {
                errormsg = "SIGMA0 can accept either one single value or as many values as the number of arguments (with NOISETYPE=MGAUSS)";
            } 
        } else {
            sigma_.resize(1, readsigma[0]);
        }

        Tools::parse(data, "SIGMA_MIN", sigma_min_);
        Tools::parse(data, "SIGMA_MAX", sigma_max_);
        Tools::parse(data, "DSIGMA", Dsigma_);
        if (Dsigma_ < 0) {
            Dsigma_ = 0.05 * (sigma_max_ - sigma_min_);
        }

        // monte carlo stuff
        Tools::parse(data, "MC_STEPS", MCsteps_);
        Tools::parse(data, "MC_STRIDE", MCstride_);

        // TODO: TESTME
        // adjust for multiple-time steps
        Tools::parse(data, "STRIDE", stride);
        MCstride_ *= stride;

        // get temperature
        double temp = 0.0;
        Tools::parse(data, "TEMP", temp);
        if (temp > 0.0) {
            kbt_ = kBoltzmann * temp;
        } else {
            kbt_ = kbt;
        }

        // variance is always the size of narg
        variance_.resize(narg, 0);

        // while sigma_mean_ has the same size of sigma
        vector<double> read_sigma_mean_;
        Tools::parseVector(data, "SIGMA_MEAN0", read_sigma_mean_);
        if (!do_optsigmamean_ && read_sigma_mean_.size() == 0 && !restart) {
            errormsg = "If you don't use OPTSIGMAMEAN and you are not RESTARTING then you MUST SET SIGMA_MEAN0";
        }

        if (noise_type_ == MGAUSS) {
            if (read_sigma_mean_.size() == narg) {
                sigma_mean_.resize(narg);
                sigma_mean_ = read_sigma_mean_;
            } else if (read_sigma_mean_.size() == 1) {
                sigma_mean_.resize(narg, read_sigma_mean_[0]);
            } else if (read_sigma_mean_.size() == 0) {
                sigma_mean_.resize(narg, 0.000001);
            } else {
                errormsg = "SIGMA_MEAN0 can accept either one single value or as many values as the arguments (with NOISETYPE=MGAUSS)";
            }
            for (unsigned i = 0; i < narg; i++) {
                if (sigma_mean_[i] > 0) {
                    variance_[i] = sigma_mean_[i] * sigma_mean_[i] * static_cast<double>(nrep_);
                }
            }
        } else {
            if (read_sigma_mean_.size() == 1) {
                sigma_mean_.resize(1, read_sigma_mean_[0]);
            } else if (read_sigma_mean_.size() == 0) {
                sigma_mean_.resize(narg, 0.000001);
            } else {
                errormsg = "If you want to use more than one SIGMA_MEAN0 you should use NOISETYPE=MGAUSS";
            }
            for (unsigned i = 0; i < narg; i++) {
                variance_[i] = sigma_mean_[0] * sigma_mean_[0] * static_cast<double>(nrep_);
            }
        } 

        // sigma mean optimisation
        if (do_optsigmamean_) {
            max_force_ = 3000.;
            Tools::parse(data, "MAX_FORCE", max_force_);
            max_force_ *= max_force_;
            sm_mod_ = 1.0;
            Tools::parse(data, "SIGMA_MEAN_MOD0", sm_mod_);
            sm_mod_min_ = 1.0;
            Tools::parse(data, "SIGMA_MEAN_MOD_MIN", sm_mod_min_);
            sm_mod_max_ = sqrt(10.);
            Tools::parse(data, "SIGMA_MEAN_MOD_MAX", sm_mod_max_);
            Dsm_mod_ = 0.01;
            Tools::parse(data, "DSIGMA_MEAN_MOD", Dsm_mod_);
        }

        if (!data.empty()) {
            errormsg = "Found the following rogue keywords in MetaInfBase input: ";
            for (unsigned i = 0; i < data.size(); ++i) {
                errormsg = errormsg + data[i] + " ";
            }
        }

        // initialize random seed
        unsigned iseed;
        if (master) {
            iseed = time(NULL)+replica_;
        } else {
            iseed = 0;     
        }
        comm.Sum(&iseed, 1);
        random[0].setSeed(-iseed);
        if (doscale_) {
            // in this case we want the same seed everywhere
            iseed = time(NULL);
            if(master) {
                multi_sim_comm.Bcast(iseed,0);
            }
            comm.Bcast(iseed, 0);
            random[1].setSeed(-iseed);
        }

        // Also seed the RNG for random index selection
        srand(time(NULL));

        // Restart from status file if necessary
        if (restart && restart_sfile.FileExist(status_file_name_)) {
            restart_sfile.open(status_file_name_);
            log.printf("  Restarting from %s\n", status_file_name_.c_str());
            double dummy;
            if (restart_sfile.scanField("time", dummy)) {
                for (unsigned i=0; i<variance_.size(); ++i) {
                    std::string msg;
                    Tools::convert(i, msg);
                    restart_sfile.scanField("variance_" + msg, variance_[i]);
                }
                for (unsigned i=0; i<sigma_.size(); ++i) {
                    std::string msg;
                    Tools::convert(i, msg);
                    restart_sfile.scanField("sigma_" + msg, sigma_[i]);
                }
                if(doscale_) {
                    restart_sfile.scanField("scale0_", scale_);
                }
                if(do_optsigmamean_) {
                    restart_sfile.scanField("sigma_mean_mod0", sm_mod_);
                }
            }
            restart_sfile.scanField();
            restart_sfile.close();
        }

        // outfile stuff
        if (write_stride_ > 0) {
            sfile.open(status_file_name_);
        }

        switch (noise_type_) {
            case GAUSS:
                log.printf("  with gaussian noise and a single noise parameter for all the data\n");
                break;
            case MGAUSS:
                log.printf("  with gaussian noise and a noise parameter for each data point\n");
                break;
            case OUTLIERS:
                log.printf("  with long tailed gaussian noise and a single noise parameter for all the data\n");
                break;
        }

        if (doscale_) {
            log.printf("  sampling a common scaling factor with:\n");
            log.printf("    initial scale parameter %f\n", scale_);
            log.printf("    minimum scale parameter %f\n", scale_min_);
            log.printf("    maximum scale parameter %f\n", scale_max_);
            log.printf("    maximum MC move of scale parameter %f\n", Dscale_);
        }

        log.printf("  number of experimental data points %u\n", narg);
        log.printf("  number of replicas %u\n", nrep_);
        log.printf("  initial data uncertainties");
        for (unsigned i=0; i<sigma_.size(); ++i) {
            log.printf(" %f", sigma_[i]);
        }
        log.printf("\n");
        log.printf("  minimum data uncertainty %f\n", sigma_min_);
        log.printf("  maximum data uncertainty %f\n", sigma_max_);
        log.printf("  maximum MC move of data uncertainty %f\n", Dsigma_);
        log.printf("  temperature of the system %f\n", kbt_);
        log.printf("  MC steps %u\n", MCsteps_);
        log.printf("  MC stride %u\n", MCstride_);
        log.printf("  initial standard errors of the mean");
        for (unsigned i=0; i<sigma_mean_.size(); ++i) {
            log.printf(" %f", sigma_mean_[i]);
        }
        log.printf("\n");
    }

    double MetaInfBase::getEnergySPE(const vector<double> &mean,
                                     const vector<double> &sigma,
                                     const double scale) {
        // calculate effective sigma
        const double smean2 = scale * scale * sigma_mean_[0] * sigma_mean_[0];
        const double s = sqrt(sigma[0] * sigma[0] + smean2);

        // cycle on arguments
        double ene = 0.0;
        for (unsigned i=0; i<narg; ++i) {
            const double dev = scale * mean[i] - parameters[i]; 
            // argument
            const double a2 = 0.5 * dev * dev + s * s;
            // increment energy
            ene += std::log(2.0 * a2 / (1.0 - exp(- a2 / smean2) ) );
        }
        // add normalization and Jeffrey's prior
        ene += std::log(sqrt(sigma[0] * sigma[0] + sigma_mean_[0] * sigma_mean_[0])) - static_cast<double>(ndata_) * std::log(sqrt2_div_pi * s);
        return kbt_ * ene;
    }

    double MetaInfBase::getEnergyGJE(const vector<double> &mean,
                                     const vector<double> &sigma,
                                     const double scale) {
        // cycle on arguments
        double ene = 0.0;
        const double scale2 = scale * scale;
        double ss = sigma[0] * sigma[0] + scale2 * sigma_mean_[0] * sigma_mean_[0];

        for (unsigned i=0; i<narg; ++i) {
            if (noise_type_ == MGAUSS) { 
                const double sigma2 = sigma[i] * sigma[i];
                const double sigma_mean2 = sigma_mean_[i] * sigma_mean_[i];
                ss = sigma2 + scale2 * sigma_mean2;
                // add Jeffrey's prior - one per sigma
                ene += 0.5 * std::log(sigma2 + sigma_mean2);
            }
            const double dev = scale * mean[i] - parameters[i];

            // deviation and normalisation 
            ene += 0.5 * dev * dev/ss + 0.5 * std::log(ss * 2. * M_PI);
        }

        // add Jeffrey's prior in case one sigma for all data points
        if (noise_type_ == GAUSS) {
            ene += 0.5 * std::log(ss);
        }
        return kbt_ * ene;
    }

    void MetaInfBase::doMonteCarlo(const vector<double> &mean_,
                                   Communicator& comm,
                                   Communicator& multi_sim_comm){

        // calculate old energy with the updated coordinates
        double old_energy = 0.;
        switch (noise_type_) {
            case GAUSS:
            case MGAUSS:
                old_energy = getEnergyGJE(mean_, sigma_, scale_);
                break;
            case OUTLIERS:
                old_energy = getEnergySPE(mean_, sigma_, scale_);
                break;
        }

        // cycle on MC steps 
        for (unsigned i=0; i<MCsteps_; ++i) {

            // propose move for scale
            double new_scale = scale_;
            if(doscale_) {
                MCtrial_++;
                if(scale_prior_ == SC_FLAT) {
                    const double r1 = random[1].Gaussian();
                    const double ds1 = sqrt(Dscale_)*r1;
                    new_scale += ds1;
                    // check boundaries
                    if(new_scale > scale_max_) {
                        new_scale = 2.0 * scale_max_ - new_scale;
                    }
                    if(new_scale < scale_min_) {
                        new_scale = 2.0 * scale_min_ - new_scale;
                    }
                } else {
                    const double r1 = random[1].Gaussian();
                    const double ds1 = (scale_mu_ - new_scale) * Dscale_ + scale_sigma_ * sqrt(2. * Dscale_) * r1;
                    new_scale += ds1;
                }

                // calculate new energy
                double new_energy = 0;
                switch (noise_type_) {
                    case GAUSS:
                    case MGAUSS:
                        new_energy = getEnergyGJE(mean_, sigma_, new_scale);
                        break;
                    case OUTLIERS:
                        new_energy = getEnergySPE(mean_, sigma_, new_scale);
                        break;
                }
                // for the scale we need to consider the total energy
                vector<double> totenergies(2);
                if (master) {
                    totenergies[0] = old_energy;
                    totenergies[1] = new_energy;
                    if (nrep_ > 1) {
                        multi_sim_comm.Sum(totenergies);
                    }
                } else {
                    totenergies[0] = 0;
                    totenergies[1] = 0;
                }
                comm.Sum(totenergies);

                // accept or reject
                const double delta = ( totenergies[1] - totenergies[0] ) / kbt_;
                // if delta is negative always accept move
                if (delta <= 0.0) {
                    old_energy = new_energy;
                    scale_ = new_scale;
                    MCaccept_++;
                    // otherwise extract random number
                } else {
                    double s = random[1].RandU01();
                    if (s < exp(-delta)) {
                        old_energy = new_energy;
                        scale_ = new_scale;
                        MCaccept_++;
                    }
                }
            }

            // propose move for sigma
            MCtrial_++;

            // Change a random element
            vector<double> new_sigma(sigma_.size());
            if (do_mc_single_) {
                new_sigma = sigma_;

                // TODO: remove slight bias towards lower values
                // choose random element
                unsigned random_index = rand() % sigma_.size();

                const double r2 = random[1].Gaussian();
                const double ds2 = sqrt(Dsigma_) * r2;
                new_sigma[random_index] = sigma_[random_index] + ds2;
                // check boundaries
                if (new_sigma[random_index] > sigma_max_) {
                    new_sigma[random_index] = 2.0 * sigma_max_ - new_sigma[random_index];
                }
                if (new_sigma[random_index] < sigma_min_) {
                    new_sigma[random_index] = 2.0 * sigma_min_ - new_sigma[random_index];
                }
            } else {
                // or change all sigmas
                for (unsigned j=0; j<sigma_.size(); j++) {
                    const double r2 = random[1].Gaussian();
                    const double ds2 = sqrt(Dsigma_) * r2;
                    new_sigma[j] = sigma_[j] + ds2;

                    // check boundaries
                    if (new_sigma[j] > sigma_max_) {
                        new_sigma[j] = 2.0 * sigma_max_ - new_sigma[j];
                    }
                    if (new_sigma[j] < sigma_min_) {
                        new_sigma[j] = 2.0 * sigma_min_ - new_sigma[j];
                    }
                }
            }

            // calculate new energy
            double new_energy=0;
            switch (noise_type_) {
                case GAUSS:
                case MGAUSS:
                    new_energy = getEnergyGJE(mean_, new_sigma, scale_);
                    break;
                case OUTLIERS:
                    new_energy = getEnergySPE(mean_, new_sigma, scale_);
                    break;
            }

            // accept or reject
            const double delta = (new_energy - old_energy) / kbt_;

            // if delta is negative always accept move
            if (delta <= 0.0) {
                old_energy = new_energy;
                sigma_ = new_sigma;
                MCaccept_++;
                // otherwise extract random number
            } else {
                const double s = random[0].RandU01();
                if (s < exp(-delta)) {
                    old_energy = new_energy;
                    sigma_ = new_sigma;
                    MCaccept_++;
                }
            }
        }

        /* save the result of the sampling */
        accept = static_cast<double>(MCaccept_) / static_cast<double>(MCtrial_);
    }

    double MetaInfBase::getEnergyForceSPE(const vector<double> &arguments,
                                          Communicator& comm,
                                          Communicator& multi_sim_comm,
                                          const vector<double> &mean,
                                          const double fact) {
        double ene = 0.0;

        const double smean2 = sigma_mean_[0] * sigma_mean_[0]; 
        const double s = sqrt(sigma_[0] * sigma_[0] + smean2);
        vector<double> f(narg, 0);

        if (master) {
            for (unsigned i=0; i<narg; ++i) {
                const double dev = scale_ * mean[i] - parameters[i]; 
                const double a2 = 0.5 * dev * dev + s * s;
                const double t = exp(-a2 / smean2);
                const double dt = 1. / t;
                const double it = 1. / (1. - t);
                const double dit = 1. / (1. - dt);
                ene += std::log(2. * a2 * it);
                f[i] = -scale_ * dev * (dit / smean2 + 1. / a2);
            }

            // collect contribution to forces and energy from other replicas
            if (nrep_ > 1) {
                multi_sim_comm.Sum(&f[0], narg);
                multi_sim_comm.Sum(&ene, 1);
            }

            // add normalizations and priors of local replica
            ene += std::log(s) - static_cast<double>(ndata_) * std::log(sqrt2_div_pi * s);
        }

        // intra-replica summation
        comm.Sum(&f[0], narg);
        comm.Sum(&ene, 1);
        double w_tmp = 0.;
        for (unsigned i=0; i<narg; ++i) {
            output_force[i] = kbt_ * fact * f[i];
            w_tmp += fact * (arguments[i] - mean[i]) * f[i];
        }

        if (do_reweight) {
            output_force[narg] = -w_tmp;
        }
        return kbt_ * ene;
    }

    double MetaInfBase::getEnergyForceGJE(const vector<double> &arguments,
                                          Communicator& comm,
                                          Communicator& multi_sim_comm,
                                          const vector<double> &mean,
                                          const double fact) {
        double ene = 0.0;

        const unsigned ssize = sigma_.size();
        vector<double> ss(ssize);
        vector<double> inv_s2(ssize, 0.);

        // if this is not MGAUSS ssize is 1
        for (unsigned i=0; i<ssize; ++i) {
            ss[i] = sigma_[i] * sigma_[i] + sigma_mean_[i] * sigma_mean_[i];
            if (master) {
                inv_s2[i] = 1.0 / ss[i];
            }
        }

        if (master && nrep_ > 1) {
            multi_sim_comm.Sum(&inv_s2[0], ssize); 
        }
        comm.Sum(&inv_s2[0],ssize);  

        double w_tmp = 0.;
        for (unsigned i=0; i<narg; ++i) {
            const double dev = scale_ * mean[i] - parameters[i]; 
            unsigned sel_sigma=0;
            if (noise_type_ == MGAUSS) {
                sel_sigma = i;
                // add Jeffrey's prior - one per sigma
                ene += 0.5 * std::log(ss[sel_sigma]);
            }
            ene += 0.5 * dev * dev * inv_s2[sel_sigma] + 0.5 * std::log(ss[sel_sigma] * 2 * M_PI);
            const double mult = fact * dev * scale_ * inv_s2[sel_sigma];
            output_force[i] = -kbt_ * mult;
            w_tmp += (arguments[i] - mean[i]) * mult;
        }

        if (do_reweight) {
            output_force[narg] = -w_tmp;
        }

        // add Jeffrey's prior in case one sigma for all data points
        if (noise_type_ == GAUSS) {
            ene += 0.5 * std::log(ss[0]);
        }
        return kbt_ * ene;
    }

    double MetaInfBase::calculate(std::vector<double>& arguments,
                                  double timestep,
                                  const long int step,
                                  const bool exchange_step,
                                  const bool checkpoint,
                                  Communicator& comm,
                                  Communicator& multi_sim_comm,
                                  OFile& sfile) {
        double norm = 0.0;
        double fact = 0.0;
        double idof = 1.0;
        double dnrep = static_cast<double>(nrep_);

        // calculate the weights either from BIAS 
        if (do_reweight) {
            vector<double> bias(nrep_, 0);
            if (master) {
                bias[replica_] = arguments[narg]; 
                if (nrep_ > 1) {
                    multi_sim_comm.Sum(&bias[0], nrep_);  
                }
            }
            comm.Sum(&bias[0], nrep_);
            const double maxbias = *(std::max_element(bias.begin(), bias.end()));
            double n2 = 0.;
            for (unsigned i=0; i<nrep_; ++i) {
                bias[i] = exp((bias[i] - maxbias) / kbt_); 
                norm += bias[i];
                n2 += bias[i] * bias[i];
            }
            fact = bias[replica_] / norm;
            idof = 1. / (1. - n2 / (norm * norm));
            // or arithmetic ones
        } else {
            norm = dnrep; 
            fact = 1.0 / norm; 
        }

        vector<double> mean(narg, 0);
        // calculate the mean 
        if (master) {
            for (unsigned i=0; i<narg; ++i) {
                mean[i] = fact * arguments[i]; 
            }
            if (nrep_ > 1) {
                multi_sim_comm.Sum(&mean[0], narg);
            }
        }
        comm.Sum(&mean[0], narg);

        if (do_optsigmamean_) {
            /* this is SIGMA_MEAN before the corrections due to the #DOF and the SCALING */
            vector<double> v_moment(narg, 0);
            if (master) {
                for (unsigned i=0; i<narg; ++i) { 
                    double tmp  = arguments[i] - mean[i];
                    v_moment[i] = fact * tmp * tmp;
                }
                if (nrep_ > 1) {
                    multi_sim_comm.Sum(&v_moment[0], narg);
                }
            }
            comm.Sum(&v_moment[0], narg);

            for (unsigned i=0; i<narg; ++i) {
                if (v_moment[i] > variance_[i]) {
                    variance_[i] = v_moment[i];
                }
            }
        }

        if (noise_type_ == MGAUSS) {
            for (unsigned i=0; i<narg; ++i) { 
                sigma_mean_[i] = sqrt(variance_[i] / dnrep);
            }
        } else {
            sigma_mean_[0] = *max_element(variance_.begin(), variance_.end());
            sigma_mean_[0] = sqrt(sigma_mean_[0] / dnrep);
        }

        /* MONTE CARLO */
        if (step % MCstride_ == 0 && !exchange_step) {
            doMonteCarlo(mean, comm, multi_sim_comm);
        }
        if (write_stride_ > 0 && (step % write_stride_ == 0 || checkpoint)) {
            writeStatus(timestep, step, sfile);
        }

        /* fix sigma_mean_ for the weighted average and the scaling factor */
        double modifier = scale_ * sqrt(idof);

        /* fix sigma_mean_ for the effect of large forces */
        if (do_optsigmamean_) {
            modifier *= sm_mod_;
        }
        for (unsigned i=0; i<sigma_mean_.size(); ++i) {
            sigma_mean_[i] *= modifier;
        }

        // calculate bias and forces
        double ene = 0; 
        switch (noise_type_) {
            case GAUSS:
            case MGAUSS:
                ene = getEnergyForceGJE(arguments, comm, multi_sim_comm, mean, fact);
                break;
            case OUTLIERS:
                ene = getEnergyForceSPE(arguments, comm, multi_sim_comm, mean, fact);
                break;
        }

        // set value of the bias
        return ene;
    }

    void MetaInfBase::writeStatus(double timestep, unsigned step, OFile& sfile_) {
        sfile_.rewind();
        sfile_.printField("time", timestep * step);
        for (unsigned i=0; i<variance_.size(); ++i) {
            std::string msg;
            Tools::convert(i,msg);
            sfile_.printField("variance_"+msg, variance_[i]);
        }
        for (unsigned i=0; i<sigma_.size(); ++i) {
            std::string msg;
            Tools::convert(i, msg);
            sfile_.printField("sigma_"+msg, sigma_[i]);
        }
        if (doscale_) {
            sfile_.printField("scale0_", scale_);
        }
        if (do_optsigmamean_) {
            sfile_.printField("sigma_mean_mod0", sm_mod_);
        }
        sfile_.printField();
        sfile_.flush();
    }

    vector<double>& MetaInfBase::getOutputForce() {
        return output_force;
    }


    MetaInfBase::MetaInfBase():
        sqrt2_div_pi(0.45015815807855),
        doscale_(false),
        scale_mu_(0),
        scale_sigma_(-1),
        scale_min_(1),
        scale_max_(-1),
        Dscale_(-1),
        Dsigma_(-1),
        sm_mod_(1.),
        ndata_(0),
        random(2),
        MCsteps_(1), 
        MCstride_(1), 
        MCaccept_(0), 
        MCtrial_(0),
        accept(0.),
        write_stride_(0),
        stride(1),
        do_reweight(false),
        do_optsigmamean_(false),
        do_mc_single_(false)
    {}

    MetaInfBase::~MetaInfBase() {
        /* if (sfile_.isOpen()) { */
        /*     sfile_.close(); */
        /* } */
    }
}
