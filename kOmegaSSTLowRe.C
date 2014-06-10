/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kOmegaSSTLowRe.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kOmegaSSTLowRe, 0);
addToRunTimeSelectionTable(RASModel, kOmegaSSTLowRe, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
tmp<volScalarField> kOmegaSSTLowRe::ReT() const
{
	tmp<volScalarField> arg
	(    	
		(k_ / (nu()*omega_))
	);
	return arg;
}



tmp<volScalarField> kOmegaSSTLowRe::alphaStar() const
{
	tmp<volScalarField> arg
	(     	
		alphaStarInf_ * ( (betaInf_ / 3.0 + (ReT()/RK_)) / ( 1.0 + (ReT()/RK_)) )
	);
	return arg;
}

tmp<volScalarField> kOmegaSSTLowRe::alpha(const volScalarField& F1) const
{
	tmp<volScalarField> arg
	(     	
		alphaInf(F1) / alphaStar() * ( (alphaZero_ + (ReT() / ROmega_)) / (1.0 + (ReT() / ROmega_)) )
	);
	return arg;
}
       
tmp<volScalarField> kOmegaSSTLowRe::betaStar() const
{
	tmp<volScalarField> arg
	(  
		betaStarInf_ *  ( (4.0/15.0 + pow4(ReT() / RBeta_)) / (1.0 + pow4(ReT() / RBeta_)) )	//non-compressible version
	);
	return arg;
}



tmp<volScalarField> kOmegaSSTLowRe::F1(const volScalarField& CDkOmega) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
        (
            max
            (
                sqrt(k_)/(0.09*omega_*y_),
                scalar(500.0)*nu()/(sqr(y_)*omega_)
            ),
            4.0*k_/(sigmaOmega2_*CDkOmegaPlus*sqr(y_))
        );

    return tanh(pow4(arg1));
}


tmp<volScalarField> kOmegaSSTLowRe::F2() const
{
    tmp<volScalarField> arg2 = max
        (
            scalar(2.0)*sqrt(k_)/(0.09*omega_*y_),
            scalar(500.0)*nu()/(sqr(y_)*omega_)
        );

    return tanh(sqr(arg2));
}


tmp<volScalarField> kOmegaSSTLowRe::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150.0*nu()/(omega_*sqr(y_)),
        scalar(10.0)
    );

    return 1.0 - tanh(pow4(arg3));
}


tmp<volScalarField> kOmegaSSTLowRe::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23() *= F3();
    }

    return f23;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kOmegaSSTLowRe::kOmegaSSTLowRe
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    betaInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaInf",
            coeffDict_,
            0.072
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            coeffDict_,
            0.0828
        )
    ),
    RBeta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "RBeta",
            coeffDict_,
            8.0
        )
    ),
    RK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "RK",
            coeffDict_,
            6.0
        )
    ),
    ROmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ROmega",
            coeffDict_,
            2.95
        )
    ),
    betaStarInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStarInf",
            coeffDict_,
            0.09
        )
    ),
    alphaStarInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaStarInf",
            coeffDict_,
            1.0
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )
    ),
    sigmaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaOmega1",
            coeffDict_,
            2.0
        )
    ),
    sigmaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaOmega2",
            coeffDict_,
            1.168
        )
    ),
    sigmaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK1",
            coeffDict_,
            1.176
        )
    ),
    sigmaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK2",
            coeffDict_,
            1.0
        )
    ),
    alphaZero_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaZero",
            coeffDict_,
            0.037
        )
    ),

    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            coeffDict_,
            false
        )
    ),

    y_(mesh_),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateOmega("omega", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    )
{
    bound(k_, kMin_);
    bound(omega_, omegaMin_);


    //not low-Re nut_, but just some wrong initialization...
    //SST:  k_/omega_ * 1/(max(1.0/alphaStar(F1),sqrt(S2)*F2()/(a1_*omega_))); cannot be used here, because F1 cannot be used
    
    nut_ =
    (
        a1_*k_
      / max
        (
            a1_*omega_,
            b1_*F2()*sqrt(2.0)*mag(symm(fvc::grad(U_)))
        )
    );
	
	Info << "------------------------------------------------------------------------" << endl;
	Info << "kOmegaSST lowRe model V1.0" << endl;
	Info << "------------------------------------------------------------------------" << endl;


  

    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> kOmegaSSTLowRe::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> kOmegaSSTLowRe::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> kOmegaSSTLowRe::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> kOmegaSSTLowRe::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


bool kOmegaSSTLowRe::read()
{
    if (RASModel::read())
    {
        beta1_.readIfPresent(coeffDict());
        beta2_.readIfPresent(coeffDict());
        a1_.readIfPresent(coeffDict());
        b1_.readIfPresent(coeffDict());
        c1_.readIfPresent(coeffDict());
        F3_.readIfPresent("F3", coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void kOmegaSSTLowRe::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
    }

    const volScalarField S2(2*magSqr(symm(fvc::grad(U_))));
    volScalarField G(GName(), nut_*S2);

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    const volScalarField CDkOmega
    (
        (2/sigmaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    const volScalarField F1(this->F1(CDkOmega));


    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::laplacian(DomegaEff(F1), omega_)
     ==
	alpha(F1)*alphaStar()*S2
      - fvm::Sp(beta(F1)*omega_, omega_)
      + fvm::SuSp
        (
            (scalar(1.0) - F1)*CDkOmega/omega_,
            omega_
        )
    );

    omegaEqn().relax();

    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omegaMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(F1), k_)
     ==
        min(G, c1_*betaStar()*k_*omega_)
      - fvm::Sp(betaStar()*omega_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity
    
    //original high Re kOmegaSST:
    //nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));
    
    // low Re kOmegaSST from paper:
    nut_ = k_/omega_ * 1.0/max(1.0/alphaStar(),sqrt(S2)*F2()/(a1_*omega_));
    
    // low Re kOmegaSST simplified:
    //nut_ = a1_*k_ /max(a1_/alphaStar(F1)*omega_,sqrt(S2)*F2());	

    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
