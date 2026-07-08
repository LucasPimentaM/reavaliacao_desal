#ifndef SOLVER

#define SOLVER

#include "../entrydata/entrydata.h"

// Defining the solver context data structure
typedef struct
{
    TS ts;
    DM da;
    EntryData entry_data;
    PetscViewer viewer;
} SolverCtx;

// Defining a solver context constructor
PetscErrorCode SolverCtxBuild(SolverCtx *solver_ctx, EntryData *entry_data);

// Defining a solver context destructor
PetscErrorCode SolverCtxDestroy(SolverCtx *solver_ctx);

#endif