#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "hac.h"


/** Hash table functions **/

/* Hash node */

typedef struct node *nodeptr;
typedef struct node {
	int dim;
	char *word;
	nodeptr left;
	nodeptr right;
} node;

/* Hash table insert */

int insert(const char *w, nodeptr *t, int d)
{
	int c;
	while ((*t) != NULL) {
		c = strcmp(w, (*t)->word);
		if (c < 0) t = &((*t)->left);
		else if (c > 0) t = &((*t)->right);
		else return -1;
	}
	*t = malloc(sizeof(node));
	(*t)->dim = d;
	(*t)->word = malloc((strlen(w)+1)*sizeof(char));
	strcpy((*t)->word, w);
	(*t)->left = (*t)->right = NULL;
	return 0;
}

/* Hash table find */

int find(const char *w, nodeptr t)
{
	int c;
	while (t != NULL) {
		c = strcmp(w, t->word);
		if (c < 0) t = t->left;
		else if (c > 0) t = t->right;
		else return t->dim;
	}
	return -1;
}

/* Hash table clear */

void del(nodeptr t)
{
	if (t == NULL) return;
	free(t->word);
	del(t->left);
	del(t->right);
	free(t);
}


/** Distance functions **/

/* Compute bytes for dimensions */

int bitcom(int d)
{
	int n, s;
	for (n = 1 ; (1<<n) != sizeof(int) ; ++n);
	s = d>>(n+3);
	return (s<<(n+3) != d ? s+1 : s);
}

/* Internal product metric */

void intprod(double **m, int **v, int l, int d)
{
	unsigned int b;
	int i, j, k, z, s = bitcom(d);
	for (i = 0 ; i != l ; ++i)
		for (j = 0 ; j <= i ; ++j) {
			for (k = z = 0 ; z != s ; ++z) {
				b = v[i][z]&v[j][z];
				while (b) {
					if (b&1) k++;
					b >>= 1;
				}
			}
			m[i][j] = m[j][i] = (double)(d-k);
		}
}

/* Jaccard distance metric */

void jaccard(double **m, int **v, int l, int d)
{
	unsigned int b;
	int i, j, k, u, z, s = bitcom(d);
	for (i = 0 ; i != l ; ++i)
		for (j = 0 ; j <= i ; ++j) {
			for (u = k = z = 0 ; z != s ; ++z) {
				b = v[i][z]&v[j][z];
				while (b) {
					if (b&1) k++;
					b >>= 1;
				}
				b = v[i][z]^v[j][z];
				while (b) {
					if (b&1) u++;
					b >>= 1;
				}
			}
			m[i][j] = m[j][i] = ((double)u)/((double)(u+k));
		}
}

/* List node */

typedef struct lnode *lnodeptr;
typedef struct lnode {
	lnodeptr next;
	int id;
} lnode;


/** Main function **/

int main(int argc, char **argv)
{
	FILE *ifp = NULL, *ofp = NULL, *cfp = NULL;
	int o, c, n, k, *r, **v;
	int i = 0, j = 0, e = 0, g = 0, z = 0, l = 0, d = 0, a = 0;
	double **m, t;
	nodeptr b = NULL;
	char *buf, *pg = argv[0];
	if (argc == 1) goto argerror;
	while (--argc)
		if (!strcmp(*(++argv), "-i")) {        /* -i for input file */
			if (!--argc || i++) goto argerror;
			ifp = fopen(*(++argv), "r");
			if (ifp == NULL) {
				if (ofp != NULL) fclose(ofp);
				if (cfp != NULL) fclose(cfp);
				fprintf(stderr, "%s: %s: Dataset file not opened\n", pg, *argv);
				return EXIT_FAILURE;
			}
		} else if (!strcmp(*argv, "-q")) {     /* -q for queue output file */
			if (!--argc || j&1) goto argerror;
			j |= 1;
			ofp = fopen(*(++argv), "w");
			if (ofp == NULL) {
				ofperror:
				if (ifp != NULL) fclose(ifp);
				if (cfp != NULL) fclose(cfp);
				fprintf(stderr, "%s: %s: Output file not created\n", pg, *argv);
				return EXIT_FAILURE;
			}
		} else if (!strcmp(*argv, "-c")) {     /* -c for cluster output file */
			if (!--argc || j&2) goto argerror;
			j |= 2;
			cfp = fopen(*(++argv), "w");
			if (cfp == NULL) goto ofperror;
		} else if (!strcmp(*argv, "-a")) {     /* -a for algorithm selection */
			if (!--argc || a) goto argerror;
			if (!strcmp(*(++argv), "s")) a = 1;
			else if (!strcmp(*argv, "c")) a = 2;
			else if (!strcmp(*argv, "a")) a = 3;
			else {
				if (ifp != NULL) fclose(ifp);
				if (ofp != NULL) fclose(ofp);
				if (cfp != NULL) fclose(cfp);
				fprintf(stderr, "%s: %s: Invalid algorithm\n", pg, *argv);
				return EXIT_FAILURE;
			}
		} else if (!strcmp(*argv, "-d")) {     /* -d for distance metric selection */
			if (!--argc || e) goto argerror;
			if (!strcmp(*(++argv), "i")) e = 1;
			else if (!strcmp(*argv, "j")) e = 2;
			else {
				if (ifp != NULL) fclose(ifp);
				if (ofp != NULL) fclose(ofp);
				if (cfp != NULL) fclose(cfp);
				fprintf(stderr, "%s: %s: Invalid metric\n", pg, *argv);
				return EXIT_FAILURE;
			}
		} else if (!strcmp(*argv, "-n")) {     /* -n for cluster number */
			if (!--argc || g) goto argerror;
			g = atoi(*(++argv));
			if (g <= 0) {
				if (ifp != NULL) fclose(ifp);
				if (ofp != NULL) fclose(ofp);
				if (cfp != NULL) fclose(cfp);
				fprintf(stderr, "%s: %s: Invalid cluster number\n", pg, *argv);
				return EXIT_FAILURE;
			}
		} else {
			argerror:
			if (ifp != NULL) fclose(ifp);
			if (ofp != NULL) fclose(ofp);
			if (cfp != NULL) fclose(cfp);
			fprintf(stderr, "%s: Invalid arguments\n", pg);
			return EXIT_FAILURE;
		}
	if (!i) {
		if (ofp != NULL) fclose(ofp);
		if (cfp != NULL) fclose(cfp);
		fprintf(stderr, "%s: No input file\n", pg);
		return EXIT_FAILURE;
	}
	if (!j) {
		if (ifp != NULL) fclose(ifp);
		fprintf(stderr, "%s: No output file\n", pg);
		return EXIT_FAILURE;
	}
	if (!a) a = 1;
	if (!e) e = 1;
	if (!g) g = 1;
	if (j&4 && g != 1) {
		if (ifp != NULL) fclose(ifp);
		if (ofp != NULL) fclose(ofp);
		if (cfp != NULL) fclose(cfp);
		fprintf(stderr, "%s: Cannot visualize\n", pg);
		return EXIT_FAILURE;
	}
	t = clock()/(double)CLOCKS_PER_SEC;            /* Read input file */
	buf = malloc(4096*sizeof(char));
	while (fgets(buf, 4096, ifp) != NULL) {
		for (i = 0 ; isspace(buf[i]) ; ++i);
		while (!isspace(buf[i])) i++;
		for (c = 0 ;; i = j, c++) {
			while (isspace(buf[i])) i++;
			if (!buf[i]) break;
			for (j = i+1 ; !isspace(buf[j]) ; ++j);
			buf[j++] = 0;
			if (!insert(&buf[i], &b, d)) d++;
		}
		if (c) l++;
	}
	for (n = 1 ; (1<<n) != sizeof(int) ; ++n);
	n += 3;
	c = d>>n;
	if (c<<n != d) c++;
	v = malloc(l*sizeof(int*));
	for (i = 0 ; i != l ; ++i) {
		v[i] = malloc(c*sizeof(int));
		for (j = 0 ; j != c ; ++j)
			v[i][j] = 0;
	}
	rewind(ifp);
	while (d && fgets(buf, 4096, ifp) != NULL) {
		for (i = 0 ; isspace(buf[i]) ; ++i);
		while (!isspace(buf[i])) i++;
		for (c = 0 ; d ; i = j, c++) {
			while (isspace(buf[i])) i++;
			if (!buf[i]) break;
			for (j = i+1 ; !isspace(buf[j]) ; ++j);
			buf[j++] = 0;
			k = find(&buf[i], b);
			if (k < 0) d = 0;
			else v[z][k>>n] |= 1<<(k&((1<<n)-1));
		}
		if (c) z++;
	}
	fclose(ifp);
	t = clock()/(double)CLOCKS_PER_SEC-t;
	free(buf);
	del(b);
	if (z != l || !d) {
		if (ifp != NULL) fclose(ifp);
		if (ofp != NULL) fclose(ofp);
		if (cfp != NULL) fclose(cfp);
		for (i = 0 ; i != l ; ++i)
			free(v[i]);
		free(v);
		fprintf(stderr, "%s: %s: Invalid file\n", argv[0], argv[1]);
		return EXIT_FAILURE;
	}
	if (t < 0.0) t = -t;
	fprintf(stderr, "Dataset file loaded           [%.2f secs]\n", t);
	fprintf(stderr, "%d vectors - %d dimensions\n", l, d);
	m = malloc(l*sizeof(double*));
	for (i = 0 ; i != l ; ++i)
		m[i] = malloc(l*sizeof(double));
	t = clock()/(double)CLOCKS_PER_SEC;            /* Compute distance matrix */
	if (e == 1) intprod(m, v, l, d);
	else jaccard(m, v, l, d);
	t = clock()/(double)CLOCKS_PER_SEC-t;
	if (t < 0.0) t = -t;
	fprintf(stderr, "Distances computed            [%.2f secs]\n", t);
	for (i = 0 ; i != l ; ++i)
		free(v[i]);
	free(v);
	r = malloc(((l-g)<<1)*sizeof(int));
	t = clock()/(double)CLOCKS_PER_SEC;
	if (a == 1) single_hac(m, r, l, g);
	else if (a == 2) complete_hac(m, r, l, g);
	else average_hac(m, r, l, g);
	t = clock()/(double)CLOCKS_PER_SEC-t;
	if (t < 0.0) t = -t;
	fprintf(stderr, "Clustering complete           [%.2f secs]\n", t);
	for (i = 0 ; i != l ; ++i)
		free(m[i]);
	free(m);
	if (ofp != NULL) {
		t = clock()/(double)CLOCKS_PER_SEC;    /* Output as connection queue */
		for (i = 0 ; i != l-g ; ++i)
			fprintf(ofp, "%d %d\n", r[i+i], r[i+i+1]);
		fclose(ofp);
		t = clock()/(double)CLOCKS_PER_SEC-t;
		fprintf(stderr, "Result queue file complete    [%.2f secs]\n", t);
	}
	if (cfp != NULL) {                             /* Output as clusters */
		lnodeptr *q, *p, u, up;
		t = clock()/(double)CLOCKS_PER_SEC;
		q = malloc(l*sizeof(lnodeptr));
		p = malloc(l*sizeof(lnodeptr));
		for (i = 0 ; i != l ; ++i) {
			q[i] = p[i] = malloc(sizeof(node));
			q[i]->next = NULL;
			q[i]->id = i;
		}
		for (i = 0 ; i != l-g ; ++i) {
			a = r[i+i];
			e = r[i+i+1];
			p[a]->next = q[e];
			p[a] = p[e];
			q[e] = NULL;
		}
		for (e = 1, i = 0 ; i != l ; ++i)
			if ((u = q[i]) != NULL) {
				fprintf(cfp, "C%d:  %d ", e++, u->id);
				for (j = 1 ; (u = u->next) != NULL ;) {
					if (!j) fputs("     ", cfp);
					fprintf(cfp, "%d ", u->id);
					if (++j == 30) {
						j = 0;
						fputc('\n', cfp);
					}
				}
				if (j) fputc('\n', cfp);
				fputc('\n', cfp);
			}
		fclose(cfp);
		t = clock()/(double)CLOCKS_PER_SEC-t;
		if (t < 0.0) t = -t;
		fprintf(stderr, "Result cluster file complete  [%.2f secs]\n", t);
		free(p);
		for (i = 0 ; i != l ; ++i)
			for (u = q[i] ; u != NULL ;) {
				up = u;
				u = u->next;
				free(up);
			}
		free(q);
	}
	free(r);
	return EXIT_SUCCESS;
}
