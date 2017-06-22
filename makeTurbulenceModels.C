// Copyright held by original authors

#include "IncompressibleTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "RASModel.H"

namespace Foam
{
    typedef IncompressibleTurbulenceModel<transportModel>
        transportModelIncompressibleTurbulenceModel;
    typedef RASModel<transportModelIncompressibleTurbulenceModel>
        RAStransportModelIncompressibleTurbulenceModel;
}

#include "kOmegaSSTLowRe.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    kOmegaSSTLowRe
);
