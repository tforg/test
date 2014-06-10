#ReactingOneDim.C documentation
####General information:

* only 1D
* **phiHsGas_**: Sensible enthalpy gas flux *J/m2/s*
* **maxDiff_**: Maximum diffussivity
* **phiGas_**: Total gas mass flux to the primary region *kg/m2/s*
* **phiHsGas_**: Sensible enthalpy gas flux *J/m2/s* [Sensible Enthalphy eqn](http://fahrzeugtechnik.fh-joanneum.at/publikationen/HAZARD_ASSESSMENT_BY_PARALLEL_CFD_MODELING/node10.html) - with negligible heat transfer through radiation, but radiation model required?!
* **chemistrySh_**: Heat release *J/s/m3*
* **Qr_**: Coupled region radiative heat flux *W/m2*. Requires user to input mapping info for coupled patches volScalarField QrCoupled_;In depth radiative heat flux *W/m2* 
*  **lostSolidMass_**: Cumulative lost mass of the condensed phase *kg*
* **addedGasMass_**: Cumulative mass generation of the gas phase *kg*
* **totalGasMassFlux_**: Total mass gas flux at the pyrolysing walls *kg/s*
* **totalHeatRR_**: Total heat release rate *J/s*
*/

'
namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{
'

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//######################################################################################################
/*
- read mesh region where pyrolysis occurs
- static
*/
//######################################################################################################

defineTypeNameAndDebug(reactingOneDim, 0);

addToRunTimeSelectionTable(pyrolysisModel, reactingOneDim, mesh);
addToRunTimeSelectionTable(pyrolysisModel, reactingOneDim, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

//######################################################################################################
/*
- read constants from system/controlDict
- read ccoeffs (radFluxName(e.g. Qr), minimumDelta) from constant/pyrolysisZones

*/
//######################################################################################################


void reactingOneDim::readReactingOneDimControls()
{
    const dictionary& solution = this->solution().subDict("SIMPLE");
    solution.lookup("nNonOrthCorr") >> nNonOrthCorr_;
    time().controlDict().lookup("maxDi") >> maxDiff_;

    coeffs().lookup("radFluxName") >> primaryRadFluxName_;
    coeffs().lookup("minimumDelta") >> minimumDelta_;
}

//######################################################################################################
/*
- check if true
*/
//######################################################################################################

}


bool reactingOneDim::read(const dictionary& dict)
{
    if (pyrolysisModel::read(dict))
    {
        readReactingOneDimControls();
        return true;
    }
    else
    {
        return false;
    }
}

//######################################################################################################
/*
- update gas flux in material from solid chemistry class
- calculate gas produced in each cell
*/
//######################################################################################################
void reactingOneDim::updatePhiGas()
{
    phiHsGas_ ==  dimensionedScalar("zero", phiHsGas_.dimensions(), 0.0);
    phiGas_ == dimensionedScalar("zero", phiGas_.dimensions(), 0.0);

    const speciesTable& gasTable = solidChemistry_->gasTable();

    forAll(gasTable, gasI)
    {
        tmp<volScalarField> tHsiGas =
            solidChemistry_->gasHs(solidThermo_.p(), solidThermo_.T(), gasI);

        const volScalarField& HsiGas = tHsiGas();

        const DimensionedField<scalar, volMesh>& RRiGas =
            solidChemistry_->RRg(gasI);

        label totalFaceId = 0;
        forAll(intCoupledPatchIDs_, i)
        {
            const label patchI = intCoupledPatchIDs_[i];

            const scalarField& phiGasp = phiHsGas_.boundaryField()[patchI];

            forAll(phiGasp, faceI)
            {
                const labelList& cells = boundaryFaceCells_[totalFaceId];
                scalar massInt = 0.0;
                forAllReverse(cells, k)
                {
                    const label cellI = cells[k];
                    massInt += RRiGas[cellI]*regionMesh().V()[cellI];
                    phiHsGas_[cellI] += massInt*HsiGas[cellI];
                }

                phiGas_.boundaryField()[patchI][faceI] += massInt;

                if (debug)
                {
                    Info<< " Gas : " << gasTable[gasI]
                        << " on patch : " << patchI
                        << " mass produced at face(local) : "
                        <<  faceI
                        << " is : " << massInt
                        << " [kg/s] " << endl;
                }
                totalFaceId ++;
            }
        }
        tHsiGas().clear();
    }
}

//######################################################################################################
/*
- update gas flux 
*/
//######################################################################################################
void reactingOneDim::updateFields()
{
    updatePhiGas();
}

//######################################################################################################
/*
- eventually move mesh 
*/
//######################################################################################################
void reactingOneDim::updateMesh(const scalarField& mass0)
{
    if (!moveMesh_)
    {
        return;
    }

    const scalarField newV(mass0/rho_);

    Info<< "Initial/final volumes = " << gSum(regionMesh().V()) << ", "
        << gSum(newV) << " [m3]" << endl;

    // move the mesh
    const labelList moveMap = moveMesh(regionMesh().V() - newV, minimumDelta_);

    // flag any cells that have not moved as non-reacting
    forAll(moveMap, i)
    {
        if (moveMap[i] == 0)
        {
            solidChemistry_->setCellReacting(i, false);
        }
    }
}

//######################################################################################################
/*
- solve Continuity after moving mesh 
-  solidChemistry_->RRg(): Return const access to the chemical source terms for gases. 
*/
//######################################################################################################
void reactingOneDim::solveContinuity()
{
    if (debug)
    {
        Info<< "reactingOneDim::solveContinuity()" << endl;
    }

    if (moveMesh_)
    {
        const scalarField mass0 = rho_*regionMesh().V();

        fvScalarMatrix rhoEqn
        (
            fvm::ddt(rho_)
         ==
          - solidChemistry_->RRg()
        );

        rhoEqn.solve();

        updateMesh(mass0);

    }
}

//######################################################################################################
/*
- solve species equation with respect to solid chemistry, moved mesh (fvc::interpolate(Yi*rho_)*regionMesh().phi())
*/
//######################################################################################################


void reactingOneDim::solveSpeciesMass()
{
    if (debug)
    {
        Info<< "reactingOneDim::solveSpeciesMass()" << endl;
    }

    volScalarField Yt(0.0*Ys_[0]);

    for (label i=0; i<Ys_.size()-1; i++)
    {
        volScalarField& Yi = Ys_[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho_, Yi)
         ==
            solidChemistry_->RRs(i)
        );

        if (regionMesh().moving())
        {
            surfaceScalarField phiYiRhoMesh
            (
                fvc::interpolate(Yi*rho_)*regionMesh().phi()
            );

            YiEqn += fvc::div(phiYiRhoMesh);

        }

        YiEqn.solve(regionMesh().solver("Yi"));
        Yi.max(0.0);
        Yt += Yi;
    }

    Ys_[Ys_.size() - 1] = 1.0 - Yt;
}

//######################################################################################################
/*
- solve energy equation (respect to moved mesh), calculate heat release + temperature of pyrolysis region
*/
//######################################################################################################

void reactingOneDim::solveEnergy()
{
    if (debug)
    {
        Info<< "reactingOneDim::solveEnergy()" << endl;
    }

    tmp<volScalarField> alpha(solidThermo_.alpha());

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho_, h_)
      - fvm::laplacian(alpha, h_)
     ==
        chemistrySh_
    );

    if (regionMesh().moving())
    {
        surfaceScalarField phihMesh
        (
            fvc::interpolate(rho_*h_)*regionMesh().phi()
        );

        hEqn += fvc::div(phihMesh);
    }

    hEqn.relax();
    hEqn.solve();

    Info<< "pyrolysis min/max(T) = " << min(solidThermo_.T()) << ", "
        << max(solidThermo_.T()) << endl;
}

//######################################################################################################
/*
- calculate mass transfer over boundaries, gas added, solid loss (integrate total heat release * time step over domain)
(- since an explicit equation can be solved immediately, the fvc namespace does so. In other words, given an operation on a volume field, the fvc namespace produces another volume field)
- calc total gas mass flux over boundaries, 
- integrate total heat release over domain
*/
//######################################################################################################

void reactingOneDim::calculateMassTransfer()
{
    totalGasMassFlux_ = 0;
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];
        totalGasMassFlux_ += gSum(phiGas_.boundaryField()[patchI]);
    }

    if (infoOutput_)
    {
        totalHeatRR_ = fvc::domainIntegrate(chemistrySh_);

        addedGasMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRg())*time_.deltaT();
        lostSolidMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRs())*time_.deltaT();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reactingOneDim::reactingOneDim(const word& modelType, const fvMesh& mesh)
:
    pyrolysisModel(modelType, mesh),
    solidChemistry_(basicSolidChemistryModel::New(regionMesh())),
    solidThermo_(solidChemistry_->solidThermo()),
    radiation_(radiation::radiationModel::New(solidThermo_.T())),
    rho_
    (
        IOobject
        (
            "rho",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidThermo_.rho()
    ),
    Ys_(solidThermo_.composition().Y()),
    h_(solidThermo_.he()),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    minimumDelta_(1e-4),

    phiGas_
    (
        IOobject
        (
            "phiGas",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    phiHsGas_
    (
        IOobject
        (
            "phiHsGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),

    chemistrySh_
    (
        IOobject
        (
            "chemistrySh",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(0.0),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0))
{
    if (active_)
    {
        read();
    }
}


reactingOneDim::reactingOneDim
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    pyrolysisModel(modelType, mesh, dict),
    solidChemistry_(basicSolidChemistryModel::New(regionMesh())),
    solidThermo_(solidChemistry_->solidThermo()),
    radiation_(radiation::radiationModel::New(solidThermo_.T())),
    rho_
    (
        IOobject
        (
            "rho",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidThermo_.rho()
    ),
    Ys_(solidThermo_.composition().Y()),
    h_(solidThermo_.he()),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    minimumDelta_(1e-4),

    phiGas_
    (
        IOobject
        (
            "phiGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    phiHsGas_
    (
        IOobject
        (
            "phiHsGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),

    chemistrySh_
    (
        IOobject
        (
            "chemistrySh",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(0.0),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0))
{
    if (active_)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

reactingOneDim::~reactingOneDim()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar reactingOneDim::addMassSources(const label patchI, const label faceI)
{
    label index = 0;
    forAll(primaryPatchIDs_, i)
    {
        if (primaryPatchIDs_[i] == patchI)
        {
            index = i;
            break;
        }
    }

    const label localPatchId =  intCoupledPatchIDs_[index];

    const scalar massAdded = phiGas_.boundaryField()[localPatchId][faceI];

    if (debug)
    {
        Info<< "\nPyrolysis region: " << type() << "added mass : "
            << massAdded << endl;
    }

    return massAdded;
}


scalar reactingOneDim::solidRegionDiffNo() const
{
    scalar DiNum = -GREAT;

    if (regionMesh().nInternalFaces() > 0)
    {
        surfaceScalarField KrhoCpbyDelta
        (
            regionMesh().surfaceInterpolation::deltaCoeffs()
          * fvc::interpolate(kappa())
          / fvc::interpolate(Cp()*rho_)
        );

        DiNum = max(KrhoCpbyDelta.internalField())*time().deltaTValue();
    }

    return DiNum;
}


scalar reactingOneDim::maxDiff() const
{
    return maxDiff_;
}


const volScalarField& reactingOneDim::rho() const
{
    return rho_;
}


const volScalarField& reactingOneDim::T() const
{
    return solidThermo_.T();
}


const tmp<volScalarField> reactingOneDim::Cp() const
{
    return solidThermo_.Cp();
}


tmp<volScalarField> reactingOneDim::kappaRad() const
{
    return radiation_->absorptionEmission().a();
}


tmp<volScalarField> reactingOneDim::kappa() const
{
    return solidThermo_.kappa();
}


const surfaceScalarField& reactingOneDim::phiGas() const
{
    return phiGas_;
}


void reactingOneDim::preEvolveRegion()
{
    pyrolysisModel::preEvolveRegion();

    // Initialise all cells as able to react
    forAll(h_, cellI)
    {
        solidChemistry_->setCellReacting(cellI, true);
    }
}


void reactingOneDim::evolveRegion()
{
    Info<< "\nEvolving pyrolysis in region: " << regionMesh().name() << endl;

    solidChemistry_->solve
    (
        time().value() - time().deltaTValue(),
        time().deltaTValue()
    );

    calculateMassTransfer();

    solveContinuity();

    chemistrySh_ = solidChemistry_->Sh()();

    updateFields();

    solveSpeciesMass();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        solveEnergy();
    }

    solidThermo_.correct();

    rho_ = solidThermo_.rho();
}


void reactingOneDim::info() const
{
    Info<< "\nPyrolysis in region: " << regionMesh().name() << endl;

    Info<< indent << "Total gas mass produced  [kg] = "
        << addedGasMass_.value() << nl
        << indent << "Total solid mass lost    [kg] = "
        << lostSolidMass_.value() << nl
        << indent << "Total pyrolysis gases  [kg/s] = "
        << totalGasMassFlux_ << nl
        << indent << "Total heat release rate [J/s] = "
        << totalHeatRR_.value() << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace regionModels
} // End namespace pyrolysisModels

// ************************************************************************* //
