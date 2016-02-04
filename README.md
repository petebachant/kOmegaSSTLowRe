kOmegaSSTLowRe turbulence model
===============================

Low Reynolds number kOmegaSST turbulence model for OpenFOAM v3.0.x.

This code was originally written by RodgriguezFatz from the cfd-online forums. See the 
original thread at 
http://www.cfd-online.com/Forums/openfoam-programming-development/134102-komegasst-lowre-damping-fluent.html


Compilation/installation
------------------------

Clone this repository into your OpenFOAM user directory and compile with `wmake`:

    cd $WM_PROJECT_USER_DIR
    git clone https://github.com/petebachant/kOmegaSSTLowRe.git
    cd kOmegaSSTLowRe
    wmake libso


Usage
-----

In `system/controlDict`:

```
libs
(
    "libmyIncompressibleRASModels.so"
);
```

In `constant/turbulenceProperties` set

    simulationType RAS;

    RAS
    {
        RASModel        kOmegaSSTLowRe;
        turbulence      on;
        printCoeffs     on;
	}


### Notes on compressibility

The current implementation only supports incompressible flow since no correction was applied to the original coded `k` and `omega` 
equations. This is due to the way `phi_` is calculated:

```c++
	const surfaceScalarField& phi_ = this->alphaRhoPhi_;
```

The model equations can be extended to a compressible formulation, although it requires rewriting the equations for `k` and `omega`, 
and the model might become invalid.


### Boundary conditions

|  Quantity | BC |
|----------:|:----|
| `p`     | `zeroGradient`  |
| `U`     |  `fixedValue (0 0 0)` | 
| `nut`   |  `nutLowReWallFunction` or `fixedValue uniform 0` |
| `k`     |  `fixedValue uniform 1e-12` |
| `omega` |  `omegaWallFunction` |
