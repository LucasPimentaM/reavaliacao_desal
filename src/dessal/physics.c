#include "../properties/properties.h"
#include "../entrydata/entrydata.h"

/*
Maxwell's model for the thermal conductivity of the membrane and empirical correlation for the Nusselt number in spacer-filled channels

Reference: I. Hitsov, K. De Sitter, C. Dotremont, P. Cauwenberg, I. Nopens, Full-scale validated Air Gap Membrane Distillation (AGMD) model
           without calibration parameters. J. Membrane Sci. 533 (2017) 309-320. https://doi.org/10.1016/j.memsci.2017.04.002
*/

PetscReal MembraneConductivity(MoistAirProperties *pore_air_prop,
                               PetscReal polymer_conductivity,
                               PetscReal membrane_porosity)
{
    PetscReal air_conductivity = pore_air_prop->thermal_conductivity, beta;

    beta = (polymer_conductivity - air_conductivity) / (polymer_conductivity + 2.0 * air_conductivity);

    return 0.93 * air_conductivity * (1.0 + 2.0 * beta * (1.0 - membrane_porosity)) / (1.0 - beta * (1.0 - membrane_porosity));
}

PetscReal ChannelHeatTransfCoef(SaltWaterProperties *bulk_water_prop,
                                SaltWaterProperties *wall_water_prop,
                                PetscReal mass_flow_rate,
                                PetscReal channel_height,
                                PetscReal channel_width,
                                PetscInt number_channels,
                                PetscReal spacer_porosity)
{
    // Properties
    PetscReal dyn_viscosity = bulk_water_prop->dyn_viscosity,
              thermal_conductivity = bulk_water_prop->thermal_conductivity,
              prandtl = bulk_water_prop->prandtl,
              wall_prandtl = wall_water_prop->prandtl;

    PetscReal mass_velocity, reynolds, nusselt;

    mass_velocity = mass_flow_rate / (number_channels * channel_height * channel_width * spacer_porosity);

    reynolds = mass_velocity * channel_height / dyn_viscosity;

    nusselt = 0.22 * PetscPowReal(reynolds, 0.69) * PetscPowReal(prandtl, 0.13);
    nusselt *= PetscPowReal(prandtl / wall_prandtl, 0.25);

    return thermal_conductivity * nusselt / channel_height;
}

PetscReal ChannelMassTransfCoef(SaltWaterProperties *bulk_water_prop,
                                SaltWaterProperties *wall_water_prop,
                                PetscReal channel_height,
                                PetscReal channel_width,
                                PetscInt number_channels,
                                PetscReal spacer_porosity)
{
    // Properties
    PetscReal dyn_viscosity = bulk_water_prop->dyn_viscosity,
              mass_diffusivity = bulk_water_prop->mass_diffusivity,
              schmidt = bulk_water_prop->schmidt,
              wall_schmidt = wall_water_prop->schmidt;

    PetscReal mass_velocity, reynolds, sherwood;

    mass_velocity = mass_flow_rate / (number_channels * channel_height * channel_width * spacer_porosity);

    reynolds = mass_velocity * channel_height / dyn_viscosity;

    sherwood = 0.22 * PetscPowReal(reynolds, 0.69) * PetscPowReal(schmidt, 0.13);
    sherwood *= PetscPowReal(schmidt / wall_schmidt, 0.25);

    return mass_diffusivity * sherwood / channel_height;
}

/*
Water mass flux across the membrane using approximate membrane and air gap permeabilities

Reference: K.M. Lisboa, D.B. Moraes, C.P. Naveira-Cotta, R.M. Cotta, Analysis of the membrane effects on the energy efficiency of water
           desalination in a direct contact membrane distillation (DCMD) system with heat recovery. Appl. Thermal Eng. 182 (2021) 116063
           https://doi.org/10.1016/j.applthermaleng.2020.116063
*/

PetscReal MolecularDiffusion(PetscReal membrane_porosity, PetscReal membrane_tortuosity, PetscReal temperature)
{
    PetscReal molecular_diffusivity;

    molecular_diffusivity = 4.46e-6 * membrane_porosity / membrane_tortuosity;
    molecular_diffusivity *= PetscPowReal(temperature, 2.334);

    return molecular_diffusivity;
}

PetscReal KnudsenDiffusion(PetscReal membrane_porosity,
                           PetscReal membrane_tortuosity,
                           PetscReal pore_diameter,
                           PetscReal temperature)
{
    PetscReal knudsen_diffusivity;

    knudsen_diffusivity = pore_diameter / 3.0;
    knudsen_diffusivity *= membrane_porosity / membrane_tortuosity;
    knudsen_diffusivity *= PetscSqrtReal(8.0 * gas_constant * temperature / (M_PI * water_molar_mass));

    return knudsen_diffusivity;
}

PetscReal SaltWaterConcentration(PetscReal mass_transfer_coef,
                                 PetscReal temperature,
                                 PetscReal salinity,
                                 PetscReal mass_flux)
{
    PetscReal molarity, density_NaCl, molar_mass_NaCl, concentration, converted_concentration;
    SaltWaterProperties prop;
    
    SaltWaterPropBuild(&prop, temperature, 0.0);

    density_NaCl = 2160.0;
    molar_mass_NaCl = 58.44e-3;
    //molarity = (1.0 / molar_mass_NaCl) * salinity / ((1.0 - salinity) / prop.density + salinity / density_NaCl); // tenho que mudar essa expressão
    /* bulk salt molarity C_b = w_b * rho / M_s  (mol/m^3) */
    molarity = salinity * prop.density / molar_mass_NaCl;
    molarity /= 1000.0; 

    concentration = molarity * PetscExpReal(mass_flux/(prop.density * mass_transfer_coef));

    converted_concentration = 1000.0 * molar_mass_NaCl * density_NaCl * concentration / (prop.density * density_NaCl + molar_mass_NaCl * concentration * (density_NaCl - prop.density));

    return converted_concentration;
}

PetscReal MassFlux(PetscReal membrane_porosity,
                   PetscReal membrane_tortuosity,
                   PetscReal membrane_thickness,
                   PetscReal pore_diameter,
                   PetscReal gap_spacer_porosity,
                   PetscReal air_gap_thickness,
                   PetscReal temperature_membrane,
                   PetscReal temperature_gap,
                   PetscReal feed_membrane_pressure,
                   PetscReal film_boundary_pressure,
                   PetscReal vacuum_pressure)
{
    PetscReal molecular_diffusivity, knudsen_diffusivity, effective_diffusivity,
              membrane_permeability, gap_permeability, permeability,
              mass_flux;
    
    temperature_membrane += 273.15;
    temperature_gap += 273.15;

    molecular_diffusivity = MolecularDiffusion(membrane_porosity, membrane_tortuosity, temperature_membrane);
    knudsen_diffusivity = KnudsenDiffusion(membrane_porosity, membrane_tortuosity, pore_diameter, temperature_membrane);

    effective_diffusivity = molecular_diffusivity * knudsen_diffusivity / (molecular_diffusivity + (atm_pressure + vacuum_pressure) * knudsen_diffusivity);

    membrane_permeability = water_molar_mass * effective_diffusivity / (gas_constant * temperature_membrane * membrane_thickness);

    molecular_diffusivity = MolecularDiffusion(1.0, 1.0, temperature_gap);

    gap_permeability = water_molar_mass * molecular_diffusivity / (gas_constant * temperature_gap * (atm_pressure + vacuum_pressure) * air_gap_thickness);

    permeability = membrane_permeability * gap_permeability / (membrane_permeability + gap_permeability);

    mass_flux = permeability * (feed_membrane_pressure - film_boundary_pressure);

    return mass_flux;
}