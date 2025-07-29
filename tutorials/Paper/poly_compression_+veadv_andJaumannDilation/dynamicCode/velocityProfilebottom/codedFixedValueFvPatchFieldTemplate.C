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
    // SHA1 = f36afe845c8c53106fb7821b305f622b42a286be
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void velocityProfilebottom_f36afe845c8c53106fb7821b305f622b42a286be(bool load)
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
    velocityProfilebottomFixedValueFvPatchVectorField
);


const char* const velocityProfilebottomFixedValueFvPatchVectorField::SHA1sum =
    "f36afe845c8c53106fb7821b305f622b42a286be";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

velocityProfilebottomFixedValueFvPatchVectorField::
velocityProfilebottomFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
{
    if (false)
    {
        Info<<"construct velocityProfilebottom sha1: f36afe845c8c53106fb7821b305f622b42a286be"
            " from patch/DimensionedField\n";
    }
}


velocityProfilebottomFixedValueFvPatchVectorField::
velocityProfilebottomFixedValueFvPatchVectorField
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
        Info<<"construct velocityProfilebottom sha1: f36afe845c8c53106fb7821b305f622b42a286be"
            " from patch/dictionary\n";
    }
}


velocityProfilebottomFixedValueFvPatchVectorField::
velocityProfilebottomFixedValueFvPatchVectorField
(
    const velocityProfilebottomFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct velocityProfilebottom sha1: f36afe845c8c53106fb7821b305f622b42a286be"
            " from patch/DimensionedField/mapper\n";
    }
}


velocityProfilebottomFixedValueFvPatchVectorField::
velocityProfilebottomFixedValueFvPatchVectorField
(
    const velocityProfilebottomFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        Info<<"construct velocityProfilebottom sha1: f36afe845c8c53106fb7821b305f622b42a286be "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

velocityProfilebottomFixedValueFvPatchVectorField::
~velocityProfilebottomFixedValueFvPatchVectorField()
{
    if (false)
    {
        Info<<"destroy velocityProfilebottom sha1: f36afe845c8c53106fb7821b305f622b42a286be\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void velocityProfilebottomFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs velocityProfilebottom sha1: f36afe845c8c53106fb7821b305f622b42a286be\n";
    }

//{{{ begin code
    #line 88 "/home/tomflint/OpenFOAM/custom/Euler_Solid_Mech/tutorials/Paper/poly_compression/0/U/boundaryField/lowerWall"
vectorField veloc(patch().size(), vector(0, 0, 0));
scalar radius=0.05e-3;
scalar Peak=-10.0;
scalar P =5.0;

            const vectorField& centre(patch().Cf());
            forAll(centre, i)
            {
                veloc[i].y() =((Peak))*Foam::exp(-Foam::pow(1.0*Foam::pow((centre[i].z()-0.1e-3),2.0)/(2.0*radius*radius),P));
            if(centre[i].z()<0.005e-3||centre[i].z()>0.195e-3)
            {
            veloc[i].y() =0.0;
            }
            }
            

operator==(veloc);
//}}} end code

    this->fixedValueFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

