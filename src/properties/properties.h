#ifndef PROPERTIES

#define PROPERTIES

#include <petscsys.h>
#include "../entrydata/entrydata.h"

// Data structure containing the thermophysical properties of moist air
typedef struct
{
    PetscReal density, specific_heat, dyn_viscosity, thermal_conductivity, prandtl;
} MoistAirProperties;

// Data structure containing the thermophysical properties of salt water
typedef struct
{
    PetscReal density, specific_heat, dyn_viscosity, thermal_conductivity, prandtl, schmidt,
              vapor_pressure, latent_heat_vaporization, mass_diffusivity;
} SaltWaterProperties;

// Function that updates the thermophysical properties of salt water
PetscErrorCode SaltWaterPropBuild(SaltWaterProperties *salt_water_prop, PetscReal temperature, PetscReal salinity);

// Function that updates the thermophysical properties of moist air
PetscErrorCode MoistAirPropBuild(MoistAirProperties *moist_air_prop, PetscReal temperature);

//Function that calculates the specific heat of dry air

PetscReal DryAirSpecificHeat(PetscReal temperature);

#endif