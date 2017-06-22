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
#include "bound.H"
#include "wallDist.H"
//#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //



// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTLowRe<BasicTurbulenceModel>::kOmegaSSTLowRe::ReT() const
{
	tmp<volScalarField> arg
	(    	
		(k_ / (this->nu()*omega_))
	);
	return arg;
}



template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTLowRe<BasicTurbulenceModel>::kOmegaSSTLowRe::alphaStar() const
{
	tmp<volScalarField> arg
	(     	
		alphaStarInf_ * ( (betaInf_ / 3.0 + (ReT()/RK_)) / ( 1.0 + (ReT()/RK_)) )
	);
	return arg;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTLowRe<BasicTurbulenceModel>::kOmegaSSTLowRe::alpha(const volScalarField& F1) const
{
	tmp<volScalarField> arg
	(     	
		alphaInf(F1) / alphaStar() * ( (alphaZero_ + (ReT() / ROmega_)) / (1.0 + (ReT() / ROmega_)) )
	);
	return arg;
}
       
template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTLowRe<BasicTurbulenceModel>::kOmegaSSTLowRe::betaStar() const
{
	tmp<volScalarField> arg
	(  
		betaStarInf_ *  ( (4.0/15.0 + pow4(ReT() / RBeta_)) / (1.0 + pow4(ReT() / RBeta_)) )	//non-compressible version
	);
	return arg;
}



template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTLowRe<BasicTurbulenceModel>::kOmegaSSTLowRe::F1(const volScalarField& CDkOmega) const
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
                scalar(500.0)*this->nu()/(sqr(y_)*omega_)
            ),
            4.0*k_/(sigmaOmega2_*CDkOmegaPlus*sqr(y_))
        );

    return tanh(pow4(arg1));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTLowRe<BasicTurbulenceModel>::kOmegaSSTLowRe::F2() const
{
    tmp<volScalarField> arg2 = max
        (
            scalar(2.0)*sqrt(k_)/(0.09*omega_*y_),
            scalar(500.0)*this->nu()/(sqr(y_)*omega_)
        );

    return tanh(sqr(arg2));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTLowRe<BasicTurbulenceModel>::kOmegaSSTLowRe::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150.0*this->nu()/(omega_*sqr(y_)),
        scalar(10.0)
    );

    return 1.0 - tanh(pow4(arg3));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTLowRe<BasicTurbulenceModel>::kOmegaSSTLowRe::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23() *= F3();
    }

    return f23;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTLowRe<BasicTurbulenceModel>::kOmegaSSTLowRe
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel> >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    betaInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaInf",
            this->coeffDict_,
            0.072
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    RBeta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "RBeta",
            this->coeffDict_,
            8.0
        )
    ),
    RK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "RK",
            this->coeffDict_,
            6.0
        )
    ),
    ROmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ROmega",
            this->coeffDict_,
            2.95
        )
    ),
    betaStarInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStarInf",
            this->coeffDict_,
            0.09
        )
    ),
    alphaStarInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaStarInf",
            this->coeffDict_,
            1.0
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    sigmaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaOmega1",
            this->coeffDict_,
            2.0
        )
    ),
    sigmaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaOmega2",
            this->coeffDict_,
            1.168
        )
    ),
    sigmaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK1",
            this->coeffDict_,
            1.176
        )
    ),
    sigmaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaZero_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaZero",
            this->coeffDict_,
            0.037
        )
    ),

    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);


    //not low-Re nut_, but just some wrong initialization...
    //SST:  k_/omega_ * 1/(max(1.0/alphaStar(F1),sqrt(S2)*F2()/(a1_*omega_))); cannot be used here, because F1 cannot be used
    
    this->nut_ =
    (
        a1_*k_
      / max
        (
            a1_*omega_,
            b1_*F2()*sqrt(2.0)*mag(symm(fvc::grad(this->U_)))
        )
    );
	
	Info << "------------------------------------------------------------------------" << endl;
	Info << "kOmegaSSTLowRe lowRe model V1.0" << endl;
	Info << "------------------------------------------------------------------------" << endl;


  

    this->nut_.correctBoundaryConditions();

    this->printCoeffs(type);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //




template<class BasicTurbulenceModel>
void kOmegaSSTLowRe<BasicTurbulenceModel>::correctNut()
{
    
}


template<class BasicTurbulenceModel>
bool kOmegaSSTLowRe<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
void kOmegaSSTLowRe<BasicTurbulenceModel>::correct()
{
    //RASModel::correct();

    if (!this->turbulence_)
    {
        return;
    }

    /*if (mesh_.changing())
    {
        y_.correct();
    }*/

    //const volScalarField S2(2*magSqr(symm(fvc::grad(this->U_))));
    const volVectorField& U = this->U_;
    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField G(this->GName(), this->nut_*S2);

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    const volScalarField CDkOmega
    (
        (2/sigmaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    const volScalarField F1(this->F1(CDkOmega));

    const surfaceScalarField& phi_ = this->alphaRhoPhi_;
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

    omegaEqn.ref().relax();

    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());

    solve(omegaEqn);
    bound(omega_, this->omegaMin_);

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

    kEqn.ref().relax();
    solve(kEqn);
    bound(k_, this->kMin_);


    // Re-calculate viscosity
    
    //original high Re kOmegaSSTLowRe:
    //nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));
    
    // low Re kOmegaSSTLowRe from paper:
    this->nut_ = k_/omega_ * 1.0/max(1.0/alphaStar(),sqrt(S2)*F2()/(a1_*omega_));
    
    // low Re kOmegaSSTLowRe simplified:
    //nut_ = a1_*k_ /max(a1_/alphaStar(F1)*omega_,sqrt(S2)*F2());	

        

    this->nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
