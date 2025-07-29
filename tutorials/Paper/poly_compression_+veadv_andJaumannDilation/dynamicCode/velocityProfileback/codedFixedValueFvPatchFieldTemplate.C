/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "codedFixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 7856814850f823b634d101147649803b95345fc3
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void velocityProfileback_7856814850f823b634d101147649803b95345fc3(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    velocityProfilebackFixedValueFvPatchVectorField
);


const char* const velocityProfilebackFixedValueFvPatchVectorField::SHA1sum =
    "7856814850f823b634d101147649803b95345fc3";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

velocityProfilebackFixedValueFvPatchVectorField::
velocityProfilebackFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
{
    if (false)
    {
        Info<<"construct velocityProfileback sha1: 7856814850f823b634d101147649803b95345fc3"
            " from patch/DimensionedField\n";
    }
}


velocityProfilebackFixedValueFvPatchVectorField::
velocityProfilebackFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct velocityProfileback sha1: 7856814850f823b634d101147649803b95345fc3"
            " from patch/dictionary\n";
    }
}


velocityProfilebackFixedValueFvPatchVectorField::
velocityProfilebackFixedValueFvPatchVectorField
(
    const velocityProfilebackFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct velocityProfileback sha1: 7856814850f823b634d101147649803b95345fc3"
            " from patch/DimensionedField/mapper\n";
    }
}


velocityProfilebackFixedValueFvPatchVectorField::
velocityProfilebackFixedValueFvPatchVectorField
(
    const velocityProfilebackFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        Info<<"construct velocityProfileback sha1: 7856814850f823b634d101147649803b95345fc3 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

velocityProfilebackFixedValueFvPatchVectorField::
~velocityProfilebackFixedValueFvPatchVectorField()
{
    if (false)
    {
        Info<<"destroy velocityProfileback sha1: 7856814850f823b634d101147649803b95345fc3\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void velocityProfilebackFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs velocityProfileback sha1: 7856814850f823b634d101147649803b95345fc3\n";
    }

//{{{ begin code
    #line 60 "/home/tomflint/OpenFOAM/custom/Euler_Solid_Mech/tutorials/Paper/poly_compression/0/U/boundaryField/Back"
vectorField veloc(patch().size(), vector(0, 0, 0));
scalar radius=0.05e-3;
scalar Peak=10.0;
scalar P =5.0;

            const vectorField& centre(patch().Cf());
            forAll(centre, i)
            {
                veloc[i].z() =((Peak))*Foam::exp(-Foam::pow(1.0*Foam::pow((centre[i].y()-0.1e-3),2.0)/(2.0*radius*radius),P));
                if(centre[i].y()<0.005e-3||centre[i].y()>0.195e-3)
            {
            veloc[i].z() =0.0;
            }
            }
            

operator==(veloc);
//}}} end code

    this->fixedValueFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

