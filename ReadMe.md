# Isca

Isca is a framework for the idealized modelling of the global circulation of
planetary atmospheres at varying levels of complexity and realism. The
framework is an outgrowth of models from GFDL designed for Earth's atmosphere,
but it may readily be extended into other planetary regimes. Various forcing
and radiation options are available. At the simple end of the spectrum a
Held-Suarez case is available. An idealized grey radiation scheme, a grey
scheme with moisture feedback, a two-band scheme and a multi-band scheme are
also available, all with simple moist effects and astronomically-based solar
forcing. At the complex end of the spectrum the framework provides a direct
connection to comprehensive atmospheric general circulation models.

For Earth modelling, options include an aqua-planet and configurable (idealized
or realistic) continents with idealized or realistic topography. Continents may
be defined by changing albedo, heat capacity and evaporative parameters, and/or
by using a simple bucket hydrology model. Oceanic Q-fluxes may be added to
reproduce specified sea-surface temperatures, with any continents or on an
aquaplanet. Planetary atmospheres may be configured by changing planetary size,
solar forcing, atmospheric mass, radiative, and other parameters.

The underlying model is written in Fortran and may largely be configured with
Python scripts, with internal coding changes required for non-standard cases.
Python scripts are also used to run the model on different architectures, to
archive the output, and for diagnostics, graphics, and post-processing. All of
these features are publicly available on a Git-based repository.

# License

Isca is distributed under a GNU GPLv3 license. See `Isca/LICENSE` file for details. 

RRTM/RRTMG: Copyright © 2002-2010, Atmospheric and Environmental Research, Inc. (AER, Inc.). 
This software may be used, copied, or redistributed as long as it is not sold and this 
copyright notice is reproduced on each copy made. This model is provided as is without 
any express or implied warranties.

The parts of ISCA provided by GFDL are also released under a GNU GPL license. A copy of the 
relevant GFDL license statment is provided below.

```
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify it    !!
!! under the terms of the GNU General Public License as published by !!
!! the Free Software Foundation, either version 3 of the License, or !!
!! (at your option) any later version.                               !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS. if not, see: http://www.gnu.org/licenses/gpl.txt  !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
```