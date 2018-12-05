[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/PetrKryslUCSD/FinEtools.jl.svg?branch=master)](https://travis-ci.org/PetrKryslUCSD/FinEtools.jl) [![codecov.io](http://codecov.io/github/PetrKryslUCSD/FinEtools.jl/coverage.svg?branch=master)](http://codecov.io/github/PetrKryslUCSD/FinEtools.jl?branch=master) 
[![Build status](https://ci.appveyor.com/api/projects/status/0qgyw2aa2529fahy?svg=true)](https://ci.appveyor.com/project/PetrKryslUCSD/finetools-jl)  [![Coverage Status](https://coveralls.io/repos/github/PetrKryslUCSD/FinEtools.jl/badge.svg?branch=master)](https://coveralls.io/github/PetrKryslUCSD/FinEtools.jl?branch=master)

# FinEtools: Finite Element tools in Julia

## News

- 11/09/2018: The name IntegData was changed to IntegDomain to better reflect the meaning of this type. Since this is an incompatible change, v1.0.0 tag was released.
- 10/02/2018: The code-coverage computation seems to be broken. The coverage in the FinEtools package hasn't actually changed and it is still at 98%.


[Past news](oldnews.md)

## Get FinEtools

This package is  registered, and hence one can do just
```julia
] add FinEtools
```
Only version 0.7, 1.0, and the nightly builds of Julia are supported. 

## Testing

```julia
] test FinEtools 
```

## Usage and Documentation

[Tutorials](https://github.com/PetrKryslUCSD/FinEtoolsTutorials.git) in the form of marked-down Julia source files using the
[Literate](https://github.com/fredrikekre/Literate.jl) workflow are available and more will  be added in the near future.

The package comes with examples  of its use 
([heat conduction](https://github.com/PetrKryslUCSD/FinEtoolsHeatConductionExamples.git), 
[linear deformation](https://github.com/PetrKryslUCSD/FinEtoolsLinearDeformationExamples.git), 
[acoustics](https://github.com/PetrKryslUCSD/FinEtoolsAcousticsExamples.git), 
[mesh generation](https://github.com/PetrKryslUCSD/FinEtoolsMeshGenerationExamples.git)). 

The documentation  is published as [Github pages](https://petrkryslucsd.github.io/FinEtools.jl). 
A use-case package explaining how to integrate FinEtools with  the user's own code is [available here](https://github.com/PetrKryslUCSD/FinEtoolsUseCase).

![Alt Visualization of sample result](http://hogwarts.ucsd.edu/~pkrysl/site.images/ScreenHunter_31%20Feb.%2009%2020.54.jpg "FinEtools.jl")
