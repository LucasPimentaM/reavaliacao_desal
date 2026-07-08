#include "../entrydata/entrydata.h"
#include "../plant/solver.h"
#include "../properties/properties.h"

PetscErrorCode ExportEntryData(EntryData *entry_data, char file[])
{
    PetscFunctionBeginUser;

    PetscViewer viewer;
    FILE *fptr;
    PetscReal membrane_area = entry_data->dessal_data.membrane_area;

    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file, &viewer);
    PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
    PetscViewerASCIIGetPointer(viewer, &fptr);

    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Entry data:,\n");
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Vacuum pressure (mbar), %.10f\n", entry_data->dessal_data.vacuum_pressure / 100.0);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Membrane area (m²), %.10f\n", membrane_area);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Membrane thickness (microns), %.10f\n", 1.0e6 * entry_data->dessal_data.membrane_thickness);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Membrane porosity (%%), %.10f\n", 100.0 * entry_data->dessal_data.membrane_porosity);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Pore diameter (microns), %.10f\n", 1.0e6 * entry_data->dessal_data.pore_diameter);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Height of the feed channel (mm), %.10f\n", 1000.0 * entry_data->dessal_data.feed_channel_height);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Height of the coolant channel (mm), %.10f\n", 1000.0 * entry_data->dessal_data.cool_channel_height);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Feed and coolant channels width (m), %.10f\n", entry_data->dessal_data.channel_width);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Number of channels, %d\n", entry_data->dessal_data.number_channels);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Channel spacer porosity (%%), %.10f\n", 100.0 * entry_data->dessal_data.spacer_porosity);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Gap spacer porosity (%%), %.10f\n", 100.0 * entry_data->dessal_data.gap_spacer_porosity);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Air gap thickness (mm), %.10f\n", 1000.0 * entry_data->dessal_data.air_gap_thickness);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Condensing wall thickness (microns), %.10f\n", 1.0e6 * entry_data->dessal_data.wall_thickness);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Thermal conductivity of the membrane material (W/mK), %.10f\n", entry_data->dessal_data.polymer_conductivity);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Thermal conductivity of the spacer material (W/mK), %.10f\n", entry_data->dessal_data.spacer_conductivity);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Dessal: Thermal conductivity of the condensing wall (W/mK), %.10f\n", entry_data->dessal_data.wall_conductivity);

    PetscViewerDestroy(&viewer);

    return 0;
}

PetscErrorCode StepMonitor(TS ts, PetscInt step, PetscReal time, Vec solution, void *ctx)
{
    SolverCtx *solver_ctx = (SolverCtx *)ctx;
    EntryData entry_data = solver_ctx->entry_data;
    DM da = solver_ctx->da;
    Vec x_local;
    PetscViewer viewer = solver_ctx->viewer;
    FILE *fptr;
    PetscScalar *x_array;
    PetscReal membrane_area = entry_data.dessal_data.membrane_area;
    PetscInt number_steps = (PetscInt)entry_data.final_time / (entry_data.num_check_steps * entry_data.init_timestep);
    SaltWaterProperties cool, feed, distilled;

    DMGetLocalVector(da, &x_local);
    DMGlobalToLocal(da, solution, INSERT_VALUES, x_local);
    DMDAVecGetArrayRead(da, x_local, &x_array);

    PetscViewerASCIIGetPointer(viewer, &fptr);

    // Salinity conversion
    SaltWaterPropBuild(&feed, x_array[0], x_array[14]);
    SaltWaterPropBuild(&cool, x_array[1], x_array[14]);
    SaltWaterPropBuild(&distilled, (x_array[4] + x_array[5]) / 2.0, 0.0);

    if (step == 0)
    {
        PetscFPrintf(PETSC_COMM_WORLD, fptr, "Time (s),"
                                             "Desal - Hot feedwater outlet temperature (°C),"
                                             "Desal - Cold feedwater outlet temperature (°C),"
                                             "Desal - Feed membrane temperature (°C),"
                                             "Desal - Membrane gap temperature (°C),"
                                             "Desal - Gap film boundary temperature (°C),"
                                             "Desal - Film wall temperature (°C),"
                                             "Desal - Coolant wall temperature (°C),"
                                             "Desal - Hot feedwater outlet salinity (wt%%),"
                                             "Desal - Hot feedwater outlet salinity (g/L),"
                                             "Desal - Cold feedwater outlet salinity (g/L),"
                                             "Desal - Mass flux (kg/m²h),"
                                             "Desal - Distilled water production rate (kg/h),"
                                             "Desal - Total heat flux (W/m²),"
                                             "Desal - Total heat transfer rate (W),"
                                             "Desal - Vapor heat flux (W/m²),"
                                             "Desal - Vapor heat transfer rate (W),"
                                             "Desal - Gain-output ratio (GOR),"
                                             "Desal - Thermal efficiency,"
                                             "Desal - Specific total energy consumption (STEC) (kWh/m^3),"
                                             "Desal - Hot feedwater outlet mass flow rate (kg/s),"
                                             "Desal - Total water produced (kg),"
                                             "Dessal - Total water produced (L),");
    }
    else if (step % number_steps == 0)
    {
        PetscFPrintf(PETSC_COMM_WORLD, fptr, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%10f,"
                                             "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",
                     time,
                     x_array[0],                                    // "Desal - Hot feedwater outlet temperature (°C),"
                     x_array[1],                                    // "Desal - Cold feedwater outlet temperature (°C),"
                     x_array[2],                                    // "Desal - Feed membrane temperature (°C),"
                     x_array[3],                                    // "Desal - Membrane gap temperature (°C),"
                     x_array[4],                                    // "Desal - Gap film boundary temperature (°C),"
                     x_array[5],                                    // "Desal - Film wall temperature (°C),"
                     x_array[6],                                    // "Desal - Coolant wall temperature (°C),"
                     100.0 * x_array[7],                            // "Desal - Hot feedwater outlet salinity (wt%%),"
                     feed.density * x_array[7],                     // "Desal - Hot feedwater outlet salinity (g/L),"
                     cool.density * x_array[14],                    // "Desal - Cold feedwater outlet salinity (g/L),"
                     3600.0 * x_array[8],                           // "Desal - Mass flux (kg/m²h),"
                     3600.0 * x_array[8] * membrane_area,           // "Desal - Distilled water production rate (kg/h),"
                     x_array[9],                                    // "Desal - Total heat flux (W/m²),"
                     x_array[9] * membrane_area,                    // "Desal - Total heat transfer rate (W),"
                     x_array[10],                                   // "Desal - Vapor heat flux (W/m²),"
                     x_array[10] * membrane_area,                   // "Desal - Vapor heat transfer rate (W),"
                     x_array[10] * membrane_area / x_array[22],     // "Desal - Gain-output ratio (GOR),"
                     100.0 * x_array[10] / x_array[9],              // "Desal - Thermal efficiency,"
                     distilled.density * x_array[53] / x_array[50], // "Desal - Specific total energy consumption (STEC) (kWh/m^3),"
                     x_array[11],                                   // "Desal - Hot feedwater outlet mass flow rate (kg/s),"
                     x_array[50],                                   // "Desal - Total water produced (kg),"
                     x_array[50] * 1000.0 / distilled.density);     // "Dessal - Total water produced (L)"
    }

    DMDAVecRestoreArrayRead(da, x_local, &x_array);
    DMRestoreLocalVector(da, &x_local);

    return 0;
}