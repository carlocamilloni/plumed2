/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#include "Colvar.h"
#include "MetaInfBase.h"
#include "ActionRegister.h"

#include <string>

using namespace std;

namespace PLMD {
    namespace colvar {

//+PLUMEDOC COLVAR CVMITest
/*
   blah
*/
//+ENDPLUMEDOC

        class CVMITest : public Colvar {
            private:
                bool pbc;
                MetaInfBase mi;
            public:
                static void registerKeywords(Keywords& keys);
                explicit CVMITest(const ActionOptions&);
                ~CVMITest();
                virtual void calculate();
        };

        PLUMED_REGISTER_ACTION(CVMITest, "CVMITEST")

        void CVMITest::registerKeywords(Keywords& keys) {
            Colvar::registerKeywords(keys);
            keys.add("atoms", "ATOMS", "involved atoms");
            keys.add("compulsory", "MI", "Metainference function");
            keys.addOutputComponent("value", "default", "42");
        }

        CVMITest::CVMITest(const ActionOptions& ao):
        PLUMED_COLVAR_INIT(ao),
        pbc(true)
        {
            bool nopbc = !pbc;
            parseFlag("NOPBC", nopbc);
            pbc = !nopbc;

            std::string midata, errors;
            parse("MI", midata);
            mi.set(midata, errors);
            if (errors.length() != 0) {
                error("Problem reading METAINFBASE keyword: " + errors);
            }

            vector<AtomNumber> atoms;
            parseAtomList("ATOMS", atoms);

            checkRead();
            addValueWithDerivatives();
            setNotPeriodic();
            requestAtoms(atoms);
        }

        CVMITest::~CVMITest() {}

        void CVMITest::calculate() {
            vector<double> positions;
            for (unsigned i = 0; i < getNumberOfAtoms(); ++i) {
                Vector pos = getPosition(i);
                positions.push_back(pos.modulo());
            }

            const double average = mi.calculate(positions);
            setBoxDerivativesNoPbc();
            setValue(average);
        }
    }
}
