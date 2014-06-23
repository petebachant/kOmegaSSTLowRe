kOmegaSSTLowRe turbulence model
===============================

Low Reynolds number kOmegaSST turbulence model for OpenFOAM.

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
    "libOpenFOAM.so"
    "libincompressibleTurbulenceModel.so"
    "libincompressibleRASModels.so"
    "libmyIncompressibleRASModels.so"
);
```

In `constant/RASProperties` set

    RASModel        kOmegaSSTLowRe;


