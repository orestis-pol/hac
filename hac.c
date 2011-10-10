#include <stdlib.h>
#include <float.h>


/*
 * Single linkage clustering
 *
 * Compute the min distance for each item
 * In each merge iteration
 *   Find the min of mins
 *   Merge the clusters
 *   Fix the min of the other clusters
 *   Fix the merged cluster
 */

void single_hac(double **a, int *r, int l, int g)
{
	int i, j, m1, m2, c, *m, *s;
	double q, **p;
	if (l <= 0 || g <= 0 || g >= l) return;
	m = malloc(l*sizeof(int));
	s = malloc(l*sizeof(int));
	p = malloc(l*sizeof(double*));
	for (i = 0 ; i != l ; ++i) {
		p[i] = malloc(l*sizeof(double));
		for (q = DBL_MAX, j = 0 ; j != l ; ++j) {
			p[i][j] = a[i][j];
			if (i != j && p[i][j] < q) {
				m[i] = j;
				q = p[i][j];
			}
		}
		s[i] = 1;
	}
	for (c = 0 ;;) {
		for (q = DBL_MAX, j = 0, i = 0 ; i != l ; ++i)
			if (m[i] >= 0 && (p[i][m[i]] < q ||
			   (p[i][m[i]] == q && s[i] < j))) {
				q = p[i][m[i]];
				j = s[i];
				m1 = i;
			}
		r[c+c] = m1;
		r[c+c+1] = m2 = m[m1];
		if (++c == l-g) break;
		m[m2] = -1;
		for (q = DBL_MAX, i = 0 ; i != l ; ++i)
			if (m[i] >= 0 && i != m1) {
				if (m[i] == m2)
					m[i] = m1;
				if (p[m2][i] < p[m1][i])
					p[m1][i] = p[m2][i];
				if (p[m1][i] < q) {
					m[m1] = i;
					q = p[m1][i];
				}
			}
		s[m1] += s[m2];
	}
	for (i = 0 ; i != l ; ++i)
		free(p[i]);
	free(p);
	free(m);
	free(s);
}


/** Auxiliary functions for heaps **/

/* Swap items */

static void swap(int *a, int i, int j)
{
	int t = a[i];
	a[i] = a[j];
	a[j] = t;
}

/* Heapify array */

static void heapify(int *h, double *v, int n)
{
	int i, c, f;
	for (h--, i = n>>1 ; i ; --i) {
		f = i;
		while ((c = (f<<1)) <= n) {
			if (c < n && v[h[c+1]] < v[h[c]]) c++;
			if (v[h[c]] >= v[h[f]]) break;
			swap(h, f, c);
			f = c;
		}
	}
}

/* Pop item from heap */

static void heap_pop(int f, int *h, double *v, int *p, int n)
{
	int c;
	swap(h, f, --n);
	swap(p, h[f], h[n]);
	while ((c = f)) {
		f = (c-1)>>1;
		if (v[h[c]] >= v[h[f]]) break;
		swap(h, f, c);
		swap(p, h[f], h[c]);
	}
	for (f = c ; (c = (f<<1)+1) < n ; f = c) {
		if (c+1 < n && v[h[c+1]] < v[h[c]]) c++;
		if (v[h[c]] >= v[h[f]]) break;
		swap(h, f, c);
		swap(p, h[f], h[c]);
	}
}

/* Push item in heap */

static void heap_push(int *h, double *v, int *p, int n)
{
	int c, f;
	for (c = n-1 ; c ; c = f) {
		f = (c-1)>>1;
		if (v[h[c]] >= v[h[f]]) break;
		swap(h, f, c);
		swap(p, h[f], h[c]);
	}
}

/*
 * Complete linkage clustering
 *
 * Create a heap for each item
 * In each merge iteration
 *   Find the min of heap tops
 *   Merge the clusters
 *   Pop the max of the two merged
 *   Fix the merged cluster
 */

void complete_hac(double **a, int *r, int l, int g)
{
	int i, j, m1, m2, c, **h, **p, *m, *s;
	double q, **v;
	if (l <= 0 || g <= 0 || g >= l) return;
	m = malloc(l*sizeof(int));
	s = malloc(l*sizeof(int));
	h = malloc(l*sizeof(int*));
	p = malloc(l*sizeof(int*));
	v = malloc(l*sizeof(double*));
	for (i = 0 ; i != l ; ++i) {
		h[i] = malloc(l*sizeof(int));
		p[i] = malloc(l*sizeof(int));
		v[i] = malloc(l*sizeof(double));
		for (j = 0 ; j != l ; ++j) {
			h[i][j] = j;
			v[i][j] = a[i][j];
		}
		heapify(h[i], v[i], l);
		for (j = 0 ; j != l ; ++j)
			p[i][h[i][j]] = j;
		heap_pop(p[i][i], h[i], v[i], p[i], l);
		m[i] = l-1;
		s[i] = 1;
	}
	for (c = 0 ;;) {
		for (q = DBL_MAX, j = l, i = 0 ; i != l ; ++i)
			if (m[i] >= 0 && (v[i][h[i][0]] < q ||
			   (v[i][h[i][0]] == q && s[i] < j))) {
				q = v[i][h[i][0]];
				j = s[i];
				m1 = i;
			}
		r[c+c] = m1;
		r[c+c+1] = m2 = h[m1][0];
		if (++c == l-g) break;
		heap_pop(0, h[m1], v[m1], p[m1], m[m1]--);
		m[m2] = -1;
		for (i = 0 ; i != l ; ++i)
			if (m[i] >= 0 && i != m1) {
				heap_pop(p[i][m2], h[i], v[i], p[i], m[i]--);
				if (v[i][m1] < v[i][m2]) {
					heap_pop(p[i][m1], h[i], v[i], p[i], m[i]);
					v[i][m1] = v[i][m2];
					heap_push(h[i], v[i], p[i], m[i]);
					heap_pop(p[m1][i], h[m1], v[m1], p[m1], m[m1]);
					v[m1][i] = v[i][m1];
					heap_push(h[m1], v[m1], p[m1], m[m1]);
				}
			}
		s[m1] += s[m2];
	}
	for (i = 0 ; i != l ; ++i) {
		free(h[i]);
		free(p[i]);
		free(v[i]);
	}
	free(h);
	free(p);
	free(v);
	free(m);
	free(s);
}

/*
 * Average linkage clustering
 *
 * Create a heap for each item
 * In each merge iteration
 *   Find the min of heap tops
 *   Merge the clusters
 *   Pop both initial distances
 *   Push new merged distance
 *   Fix the merged cluster
 */

void average_hac(double **a, int *r, int l, int g)
{
	int i, j, m1, m2, c, **h, **p, *m, *s;
	double q, **v;
	if (l <= 0 || g <= 0 || g >= l) return;
	m = malloc(l*sizeof(int));
	s = malloc(l*sizeof(int));
	h = malloc(l*sizeof(int*));
	p = malloc(l*sizeof(int*));
	v = malloc(l*sizeof(double*));
	for (i = 0 ; i != l ; ++i) {
		h[i] = malloc(l*sizeof(int));
		p[i] = malloc(l*sizeof(int));
		v[i] = malloc(l*sizeof(double));
		for (j = 0 ; j != l ; ++j) {
			h[i][j] = j;
			v[i][j] = a[i][j];
		}
		heapify(h[i], v[i], l);
		for (j = 0 ; j != l ; ++j)
			p[i][h[i][j]] = j;
		heap_pop(p[i][i], h[i], v[i], p[i], l);
		m[i] = l-1;
		s[i] = 1;
	}
	for (c = 0 ;;) {
		for (q = DBL_MAX, j = l, i = 0 ; i != l ; ++i)
			if (m[i] >= 0 && (v[i][h[i][0]] < q ||
			   (v[i][h[i][0]] == q && s[i] < j))) {
				q = v[i][h[i][0]];
				j = s[i];
				m1 = i;
			}
		r[c+c] = m1;
		r[c+c+1] = m2 = h[m1][0];
		if (++c == l-g) break;
		heap_pop(0, h[m1], v[m1], p[m1], m[m1]--);
		m[m2] = -1;
		for (i = 0 ; i != l ; ++i)
			if (m[i] >= 0 && i != m1) {
				heap_pop(p[i][m2], h[i], v[i], p[i], m[i]--);
				heap_pop(p[i][m1], h[i], v[i], p[i], m[i]);
				v[i][m1] = (v[i][m1]*(double)s[m1]+v[i][m2]*(double)s[m2])/
				           (double)(s[m1]+s[m2]);
				heap_push(h[i], v[i], p[i], m[i]);
				heap_pop(p[m1][i], h[m1], v[m1], p[m1], m[m1]);
				v[m1][i] = v[i][m1];
				heap_push(h[m1], v[m1], p[m1], m[m1]);
			}
		s[m1] += s[m2];
	}
	for (i = 0 ; i != l ; ++i) {
		free(h[i]);
		free(p[i]);
		free(v[i]);
	}
	free(h);
	free(p);
	free(v);
	free(m);
	free(s);
}
