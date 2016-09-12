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
        keys.add("compulsory", "TEST", "testing");
    }

    void MetaInfBase::set(const std::string& definition, std::string& errormsg) {
        vector<string> data = Tools::getWords(definition);
        if (data.size() < 1) {
            errormsg = "Missing all input for MetaInfBase";
        }

        test = 0.0;
        Tools::parse(data, "TEST", test);

        if (!data.empty()) {
            errormsg = "Found the following rogue keywords in MetaInfBase input: ";
            for (unsigned i = 0; i < data.size(); ++i) {
                errormsg = errormsg + data[i] + " ";
            }
        }
    }

    double MetaInfBase::calculate(std::vector<double>& arguments) {
        const unsigned ndata = arguments.size();
        double tot = 0.0;
        for (unsigned i = 0; i < ndata; ++i) {
            tot += arguments[i];
        }
        return tot / static_cast<double>(ndata);
    }

    MetaInfBase::MetaInfBase():
        test(0.0)
    {}

    MetaInfBase::~MetaInfBase() {}
}
