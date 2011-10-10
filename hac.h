#ifndef HAC
#define HAC


/*
 * Single complete and average linkage
 * hierarchical agglomerative clustering
 *
 *  First argument is the table of distances
 *  Second argument is the result merge queue
 *  Third argument is the number of items
 *  Fourth argument is the number of clusters
 */

void single_hac(double **p, int *r, int l, int g);

void complete_hac(double **p, int *r, int l, int g);

void average_hac(double **p, int *r, int l, int g);


#endif
