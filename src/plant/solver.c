#include "solver.h"

#define NUM_VAR 14

PetscErrorCode SolverCtxBuild(SolverCtx *solver_ctx, EntryData *entry_data)
{
    PetscFunctionBeginUser;
    TS ts; // passos de tempo
    TSAdapt adapt;
    SNES snes; // não linear solver
    KSP ksp;   // linear solver
    DM da;     // forma a malha
    PetscViewer viewer;

    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, NUM_VAR, 1, 1, NULL, &da);
    DMSetUp(da);

    TSCreate(PETSC_COMM_WORLD, &ts);
    TSSetProblemType(ts, TS_NONLINEAR);
    TSSetEquationType(ts, TS_EQ_DAE_SEMI_EXPLICIT_INDEX1);
    TSSetType(ts, TSBEULER);
    TSSetDM(ts, da);
    TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP);
    TSSetTimeStep(ts, entry_data->init_timestep);
    TSSetMaxTime(ts, entry_data->final_time);
    TSGetAdapt(ts, &adapt);
    TSAdaptSetType(adapt, TSADAPTNONE);
    TSGetSNES(ts, &snes);

    SNESSetType(snes, SNESNEWTONLS);
    SNESSetTolerances(snes, 1.0e-8, 1.0e-8, PETSC_DEFAULT, 10000, -1);
    SNESGetKSP(snes, &ksp);
    KSPSetTolerances(ksp, 1.0e-8, 1.0e-8, PETSC_DEFAULT, 10000);
    KSPGMRESSetOrthogonalization(ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);

    TSSetFromOptions(ts);

    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "./results/report.csv", &viewer);
    PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);

    solver_ctx->ts = ts;
    solver_ctx->da = da;
    solver_ctx->entry_data = *entry_data;
    solver_ctx->viewer = viewer;

    return 0;
}

PetscErrorCode SolverCtxDestroy(SolverCtx *solver_ctx)
{
    PetscFunctionBeginUser;
    TSDestroy(&solver_ctx->ts);
    DMDestroy(&solver_ctx->da);
    PetscViewerDestroy(&solver_ctx->viewer);

    return 0;
}