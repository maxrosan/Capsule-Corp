#include "mymath.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

void v3d_set(v3d *v, double x, double y, double z) {

#ifdef DEBUG
	assert(v != NULL);
#endif

	v->x = x;
	v->y = y;
	v->z = z;
}

void v3d_minus(v3d *a, v3d *b, v3d *res) {

#ifdef DEBUG
	assert(a != NULL);
	assert(b != NULL);
	assert(res != NULL);
#endif

	res->x = a->x - b->x;	res->y = a->y - b->y;
	res->z = a->z - b->z;
}


double v3d_length(v3d *v) {

#ifdef DEBUG
	assert(v != NULL);
#endif

	return sqrt(
	 (v->x * v->x) +
	 (v->y * v->y) +
	 (v->z * v->z));
}

void v3d_copy(v3d *res, v3d *a) {

#ifdef DEBUG	
	assert(res != NULL);
	assert(a != NULL);
#endif

	v3d_set(res, a->x, a->y, a->z);
}

void v3d_cross(v3d *a, v3d *b, v3d *res) {

#ifdef DEBUG
	assert(a != NULL);
	assert(b != NULL);
	assert(res != NULL);
#endif

	res->x = a->y * b->z - a->z * b->y;
	res->y = a->z * b->x - a->x * b->z;
	res->z = a->x * b->y - a->y * b->x;
}

void v3d_print(v3d *v) {

#ifdef DEBUG
	assert(v != NULL);
#endif

	printf("(%lf, %lf, %lf)", v->x, v->y, v->z);
}

double v3d_dot(v3d *a, v3d *b) {

#ifdef DEBUG
	assert(a != NULL);
	assert(b != NULL);
#endif

	return (
	 a->x * b->x +
	 a->y * b->y +
	 a->z * b->z);
}

void v3d_normalize(v3d *a) {
	
	double mod;

#ifdef DEBUG
	assert(a != NULL);
#endif
	
	mod = v3d_length(a);

#ifdef DEBUG
	assert(fabs(mod) >= 1e-10);
#endif
	
	a->x /= mod;
	a->y /= mod;
	a->z /= mod;
}
