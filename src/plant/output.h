#ifndef OUTPUT

#define OUTPUT

#include "../entrydata/entrydata.h"

// Function to export the entry data to a file
PetscErrorCode ExportEntryData(EntryData *entry_data, char file[]);

// Function to monitor the steps into a timeseries file
PetscErrorCode StepMonitor(TS ts, PetscInt step, PetscReal time, Vec solution, void *ctx);

#endif