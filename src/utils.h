#ifndef R2D2_UTILS_INCLUDED
#define R2D2_UTILS_INCLUDED

void rank_bycol(int nrow, int ncol, double *datamat, int*rankmat);

void sort_bycol(int nrow, int ncol, double* datamat, double* sortedmat);

void scale_rank(int nrow, int ncol, int* mat, double scale_factor);

void shrink_bcsorted(int nrow_ref, int nrow_bc, int ncol_bc, double *bc_sorted[]);

void find_closestrank(
    int* bestmatch,
    double* bestdistance,
    int isearchstart,
    int ntimesteps,
    int nconddim,
    int iconddim[],
    int nrow_bc,
    int ranks_bc[],
    int nrow_ref,
    int ncol_ref,
    int ranks_ref[]
);

void search_byblock(
    int* nsearch,
    double  r2d2_bc[],
    int *time_bestanalogue[],
    double *dist_bestanalogue[],
    int visited_time[],
    int nconddim,
    int iconddim[],
    int nrow_bc,
    int ranks_bc[],
    double sorted_bc[],
    int nrow_ref,
    int ncol_ref,
    int ranks_ref[],
    int lag_search,
    int lag_keep
);

#endif
