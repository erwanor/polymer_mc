#include "tensor3.h"

dtensor3 dtensor3_dot(const dvec* v0,
	const dvec* v1)
{
	dtensor3 ret;
	ret.xx = v0->x * v1->x; ret.xy = v0->x * v1->y; ret.xz = v0->x * v1->z;
	ret.yx = v0->y * v1->x; ret.yy = v0->y * v1->y; ret.yz = v0->y * v1->z;
	ret.zx = v0->z * v1->x; ret.zy = v0->z * v1->y; ret.zz = v0->z * v1->z;
	return ret;
}

void dtensor3_add(dtensor3* t0,
	const dtensor3* t1)
{
	t0->xx += t1->xx; t0->xy += t1->xy; t0->xz += t1->xz;
	t0->yx += t1->yx; t0->yy += t1->yy; t0->yz += t1->yz;
	t0->zx += t1->zx; t0->zy += t1->zy; t0->zz += t1->zz;
}

dtensor3 dtensor3_add_new(const dtensor3* t0,
	const dtensor3* t1)
{
	dtensor3 ret;
	ret.xx = t0->xx + t1->xx; ret.xy = t0->xy + t1->xy; ret.xz = t0->xz + t1->xz;
	ret.yx = t0->yx + t1->yx; ret.yy = t0->yy + t1->yy; ret.yz = t0->yz + t1->yz;
	ret.zx = t0->zx + t1->zx; ret.zy = t0->zy + t1->zy; ret.zz = t0->zz + t1->zz;
	return ret;
}

void dtensor3_sub(dtensor3* t0,
	const dtensor3* t1)
{
	t0->xx -= t1->xx; t0->xy -= t1->xy; t0->xz -= t1->xz;
	t0->yx -= t1->yx; t0->yy -= t1->yy; t0->yz -= t1->yz;
	t0->zx -= t1->zx; t0->zy -= t1->zy; t0->zz -= t1->zz;
}

dtensor3 dtensor3_sub_new(const dtensor3* t0,
	const dtensor3* t1)
{
	dtensor3 ret;
	ret.xx = t0->xx - t1->xx; ret.xy = t0->xy - t1->xy; ret.xz = t0->xz - t1->xz;
	ret.yx = t0->yx - t1->yx; ret.yy = t0->yy - t1->yy; ret.yz = t0->yz - t1->yz;
	ret.zx = t0->zx - t1->zx; ret.zy = t0->zy - t1->zy; ret.zz = t0->zz - t1->zz;
	return ret;
}