# crystaldynamicsFoam

An open-source, cell-centred finite volume solver for the solid dynamics of
polycrystalline, multi-material solids under extreme loading. The solver is
fully Eulerian: the grid stays fixed, the material and all of its history
(stresses, lattice orientations, hardening state, microstructure) is advected
through it. This mean there is no mesh distortion problem when the deformation
localise, which is the usual failure mode of Lagrangian crystal plasticity
codes.



## Authors
#### Tom Flint (University of Manchester)
#### Pratheek Shanthraj (United Kingdom Atomic Energy Authority)
#### Philip Cardiff (University College Dublin)
#### Joao Quinta da Fonseca (University of Manchester)
#### Joe Robson (University of Manchester)
#### Harry Smith (University of Manchester)



## Installation

The solver is built on top of the OpenFOAM libraries (openfoam.com, tested
from v2506 onwards) and uses the header-only
[Fastor](https://github.com/romeric/Fastor) library for the rank-four tensor
algebra. It is the first solver of a more general solid-dynamics toolbox,
updated-Lagrangian and material-point-method solvers will follow, sharing the
same material and interface infrastructure.

We like the ROUND advection schemes available at [libROUNDSchemes](https://github.com/ROUNDschemes/libROUNDSchemes). Most of the tutorials use these schemes so you should download this library and run the top-level allwmake command after sourcing the above referenced openfoam environment.

The full description of the model and the numerics is in the papers (see
[Citing](#citing--licence) below), here we give only the physics at top level.

## Physics

**Conservation laws.** Mass, momentum and energy are solved in conservative
form on the fixed grid,

$$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \boldsymbol{U}) = 0,
$$
 
$$
\frac{\partial (\rho \boldsymbol{U})}{\partial t} +
\nabla \cdot (\rho \boldsymbol{U} \otimes \boldsymbol{U}) -
\nabla \cdot \boldsymbol{\sigma} = 0,
$$

$$
\frac{\partial (\rho E)}{\partial t} +
\nabla \cdot (\rho \boldsymbol{U} E) -
\nabla \cdot (\boldsymbol{\sigma} \cdot \boldsymbol{U}) -
\nabla \cdot (\kappa \nabla T) = 0,
$$

with $E = e + \boldsymbol{U}^2/2$. All history variables are transported in
the conservative form also, using high-resolution (ROUND) convection schemes
to keep the advected profiles sharp.

**Crystal plasticity.** Each phase $i$ evolve its own Cauchy stress through
the Jaumann rate with respect to the lattice spin,

$$
\dot{\boldsymbol{\sigma}}_i =
\mathbb{C}_i : (\boldsymbol{D} - \boldsymbol{D}^p_i)
- \boldsymbol{\sigma}_i \, \mathrm{tr}(\boldsymbol{D} - \boldsymbol{D}^p_i)
+ (\boldsymbol{\Omega} - \boldsymbol{\Omega}^p_i) \cdot \boldsymbol{\sigma}_i
- \boldsymbol{\sigma}_i \cdot (\boldsymbol{\Omega} - \boldsymbol{\Omega}^p_i),
$$

where $\mathbb{C}_i$ is the anisotropic stiffness rotated to the current
lattice orientation (stored as Euler angles and updated by an exponential map
of the elastic spin). Plastic flow is by slip on the user-defined slip
systems, with the rate-dependent power law

$$
\dot{\gamma}^{\alpha}_i = \dot{\gamma}_0
\left| \frac{\tau^{\alpha}_i}{g^{\alpha}_i} \right|^{k}
\mathrm{sgn}(\tau^{\alpha}_i),
\qquad \tau^{\alpha}_i = \boldsymbol{\sigma} : \boldsymbol{p}^{\alpha}_i,
$$

and a Voce-type saturating hardening law for the slip resistances
$g^{\alpha}_i$. The mixture stress that enters the momentum equation is the
order-parameter weighted sum, $\boldsymbol{\sigma} = \sum_i \varphi_i
\boldsymbol{\sigma}_i$.

**Microstructure.** Grains and phases are represented by diffuse-interface
order parameters $\varphi_i \in [0,1]$, $\sum_i \varphi_i = 1$. They evolve
with a multi-phase-field equation driven by the interfacial energy and by the
stored elastic energy of each phase,

$$
\dot{\varphi}_i = -\frac{1}{\tilde{N}} \sum_{j \neq i} M_{ij}
\left( \frac{\delta F}{\delta \varphi_i} - \frac{\delta F}{\delta \varphi_j} \right),
\qquad
F = F_{\mathrm{surf}} + F_{\mathrm{mech}},
$$

where the mechanical driving force per phase is the elastic energy density
$\tfrac{1}{2} \boldsymbol{\sigma}_i : \mathbb{C}_i^{-1} :
\boldsymbol{\sigma}_i$. So the stress state drives the microstructure
evolution directly (interface migration, growth of damage). Setting the
mobilities to zero, or `PF_drivingForceActive false;`, recovers the pure advection of the
microstructure.

**Shock compression (optional).** For the strong shocks, where the
linear-elastic volumetric response is not valid anymore, the solver provides a
Mie–Grüneisen equation of state (`EOSMode MieGruneisen;`). The stress is then
split: crystal plasticity evolves only the per-phase deviators and the
pressure comes algebraically from the conserved density and internal energy,

$$
\boldsymbol{\sigma} = \sum_i \varphi_i \boldsymbol{\sigma}^{d}_i - (p + q)\boldsymbol{I},
\qquad
p(\rho, e) =
\begin{cases}
p_H(\mu) + \Gamma_0 \rho_0 \left[ e - e_H(\mu) \right], & \mu \ge 0,\\
\rho_0 c_0^2 \mu + \Gamma_0 \rho_0 \left[ e - \tfrac{1}{2} c_0^2 \mu^2 \right], & \mu < 0,
\end{cases}
$$

with the compression $\mu = \rho/\rho_0 - 1$ and $p_H$, $e_H$ the principal
Hugoniot implied from $U_s = c_0 + s\,u_p$. Shocks are captured with the
standard von Neumann–Richtmyer artificial viscosity $q$, and the shock
heating feeds back to the pressure through the energy equation, no extra
source terms are needed.

## Inputs

All the physics is configured from a single dictionary,
`constant/GrainProperties`, no code changes are needed for different crystal
structures or numbers of phases. The top-level entries:

| Keyword | Type | Meaning |
|---|---|---|
| `N_Phases` | label | Number of order parameters |
| `N_Materials` | label | Number of material definitions |
| `Phase2Material` | labelList | Material index of each order parameter |
| `n_constitutive_substeps` | label | Sub-steps of the constitutive integration |
| `phi1`, `phi2`, `phi3` | scalarList | Initial Euler angles per order parameter (rad) |
| `interface_width` | scalar | Diffuse interface width (m) |
| `Energy`, `Mobility` | hashTable | Pairwise interfacial energies and mobilities |
| `stressUpdate` | Switch | Optional; `false` gives the prescribed-stress mode (default `true`) |
| `PerPhaseSigma` | Switch | Optional; distinct constant stress per order parameter in the prescribed-stress mode (default `false`) |
| `sigmaPrescribed` | List\<tensor\> | Required if `PerPhaseSigma true;`, one Cauchy stress per order parameter |
| `PF_drivingForceActive` | Switch | Optional; activates the phase-field driving force (default `false`) |
| `verbose` | Switch | Optional; tensor-inversion self-tests and exit (default `false`) |
| `EOSMode` | word | Optional; `linearElastic` (default) or `MieGruneisen` |
| `EOS_Tref` | scalar | Optional; EOS reference temperature (default 300 K) |
| `EOS_qQuad`, `EOS_qLinear` | scalar | Optional; artificial viscosity coefficients (defaults 1.2, 0.06) |
| `EOS_waveStabiliser` | scalar | Optional; momentum stabiliser coefficient (default 0, off) |
| `PF_maxDeltaPhi` | scalar | Optional; per-step order-parameter change cap with `PF_drivingForceActive true;` (default 0.05) |

And per material, in the sub-dictionaries `Material.<i>`:

| Keyword | Type | Meaning |
|---|---|---|
| `Cmnpq` | scalar | Stiffness components in the crystal frame, zero-based indices (Pa) |
| `rho_0` | scalar | Reference density (kg m⁻³) |
| `Cv` | scalar | Specific heat capacity (J kg⁻¹ K⁻¹) |
| `kappa` | scalar | Thermal conductivity (W m⁻¹ K⁻¹) |
| `N_slips` | label | Number of slip systems |
| `m<j>`, `n<j>` | vector | Slip direction and plane normal of system j (crystal frame) |
| `initialCRSS` | scalar | Initial CRSS (Pa) |
| `tau_S` | scalar | Saturation shear stress (Pa) |
| `h_0` | scalar | Initial hardening rate (Pa) |
| `r` | scalar | Latent hardening parameter |
| `k_exponent` | scalar | Rate-sensitivity exponent |
| `gammadot_0` | scalar | Reference slip rate (s⁻¹) |
| `c0` | scalar | Bulk sound speed (m s⁻¹); required for `MieGruneisen` |
| `s_UsUp` | scalar | Hugoniot slope; required for `MieGruneisen` |
| `Gamma0` | scalar | Grüneisen parameter; required for `MieGruneisen` |

## Building and running

```
./Allwmake
```

Needs a recent OpenFOAM (openfoam.com, v2506 onwards) sourced, and the Fastor
include path set in `Make/options`. The ROUND convection schemes library is
optional but recommended, it improves a lot the advection accuracy of the
history fields. Cases run like any OpenFOAM solver (`blockMesh`, set the
fields, `crystaldynamicsFoam`; for parallel `decomposePar` and `mpirun
-np <n> crystaldynamicsFoam -parallel`).

The `tutorials/` folder contains the verification cases from the paper
(elastic waves, planar-shock Hugoniot, Riemann problems, phase-field interface
relaxation) plus the demonstration cases.

Some of these tutorials contain a `dynamicCode` directory,
and that code also neede to be compiled, using `wmake libso .` from within each
subdirectory e.g. `tutorials/crystaldynamicsFoam/Circle_Squash/dynamicCode/velocityProfileback`.

## Citing / Licence

GPLv3. If you use the code please cite the methodology paper, T. Flint et
al., *International Journal of Mechanical Sciences* 325 (2026) 111794, and
the software paper (submitted to *Computer Physics Communications*).
