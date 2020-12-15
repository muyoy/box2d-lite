/*
* Copyright (c) 2006-2007 Erin Catto http://www.gphysics.com
*
* Permission to use, copy, modify, distribute and sell this software
* and its documentation for any purpose is hereby granted without fee,
* provided that the above copyright notice appear in all copies.
* Erin Catto makes no representations about the suitability 
* of this software for any purpose.  
* It is provided "as is" without express or implied warranty.
*/

#include "box2d-lite/Arbiter.h"
#include "box2d-lite/Body.h"
#include <stdio.h>

// Box vertex and edge numbering:
//
//        ^ y
//        |
//        e1
//   v2 ------ v1
//    |        |
// e2 |        | e4  --> x
//    |        |
//   v3 ------ v4
//        e3

enum Axis
{
	FACE_A_X,
	FACE_A_Y,
	FACE_B_X,
	FACE_B_Y
};

enum EdgeNumbers
{
	NO_EDGE = 0,
	EDGE1,
	EDGE2,
	EDGE3,
	EDGE4
};

struct ClipVertex
{
	ClipVertex() { fp.value = 0; }
	Vec2 v;
	FeaturePair fp;
};

void Flip(FeaturePair& fp)
{
	Swap(fp.e.inEdge1, fp.e.inEdge2);
	Swap(fp.e.outEdge1, fp.e.outEdge2);
}

int ClipSegmentToLine(ClipVertex vOut[2], ClipVertex vIn[2],
					  const Vec2& normal, float offset, char clipEdge)
{
	// Start with no output points
	int numOut = 0;

	// Calculate the distance of end points to the line
	float distance0 = Dot(normal, vIn[0].v) - offset;
	float distance1 = Dot(normal, vIn[1].v) - offset;

	// If the points are behind the plane
	if (distance0 <= 0.0f) vOut[numOut++] = vIn[0];
	if (distance1 <= 0.0f) vOut[numOut++] = vIn[1];

	// If the points are on different sides of the plane
	if (distance0 * distance1 < 0.0f)
	{
		// Find intersection point of edge and plane
		float interp = distance0 / (distance0 - distance1);
		vOut[numOut].v = vIn[0].v + interp * (vIn[1].v - vIn[0].v);
		if (distance0 > 0.0f)
		{
			vOut[numOut].fp = vIn[0].fp;
			vOut[numOut].fp.e.inEdge1 = clipEdge;
			vOut[numOut].fp.e.inEdge2 = NO_EDGE;
		}
		else
		{
			vOut[numOut].fp = vIn[1].fp;
			vOut[numOut].fp.e.outEdge1 = clipEdge;
			vOut[numOut].fp.e.outEdge2 = NO_EDGE;
		}
		++numOut;
	}

	return numOut;
}

static void ComputeIncidentEdge(ClipVertex c[2], const Vec2& h, const Vec2& pos,
								const Mat22& Rot, const Vec2& normal)
{
	// The normal is from the reference box. Convert it
	// to the incident boxe's frame and flip sign.
	Mat22 RotT = Rot.Transpose();
	Vec2 n = -(RotT * normal);
	Vec2 nAbs = Abs(n);

	if (nAbs.x > nAbs.y)
	{
		if (Sign(n.x) > 0.0f)
		{
			c[0].v.Set(h.x, -h.y);
			c[0].fp.e.inEdge2 = EDGE3;
			c[0].fp.e.outEdge2 = EDGE4;

			c[1].v.Set(h.x, h.y);
			c[1].fp.e.inEdge2 = EDGE4;
			c[1].fp.e.outEdge2 = EDGE1;
		}
		else
		{
			c[0].v.Set(-h.x, h.y);
			c[0].fp.e.inEdge2 = EDGE1;
			c[0].fp.e.outEdge2 = EDGE2;

			c[1].v.Set(-h.x, -h.y);
			c[1].fp.e.inEdge2 = EDGE2;
			c[1].fp.e.outEdge2 = EDGE3;
		}
	}
	else
	{
		if (Sign(n.y) > 0.0f)
		{
			c[0].v.Set(h.x, h.y);
			c[0].fp.e.inEdge2 = EDGE4;
			c[0].fp.e.outEdge2 = EDGE1;

			c[1].v.Set(-h.x, h.y);
			c[1].fp.e.inEdge2 = EDGE1;
			c[1].fp.e.outEdge2 = EDGE2;
		}
		else
		{
			c[0].v.Set(-h.x, -h.y);
			c[0].fp.e.inEdge2 = EDGE2;
			c[0].fp.e.outEdge2 = EDGE3;

			c[1].v.Set(h.x, -h.y);
			c[1].fp.e.inEdge2 = EDGE3;
			c[1].fp.e.outEdge2 = EDGE4;
		}
	}

	c[0].v = pos + Rot * c[0].v;
	c[1].v = pos + Rot * c[1].v;
}
int BNBCollide(Contact*contacts, Body* bodyA, Body* bodyB)
{
	// Setup
	Vec2 hA = 0.5f * bodyA->width;
	Vec2 hB = 0.5f * bodyB->width;

	Vec2 posA = bodyA->position;
	Vec2 posB = bodyB->position;

	Mat22 RotA(bodyA->rotation), RotB(bodyB->rotation);

	Mat22 RotAT = RotA.Transpose();
	Mat22 RotBT = RotB.Transpose();

	Vec2 dp = posB - posA;
	Vec2 dA = RotAT * dp;
	Vec2 dB = RotBT * dp;

	Mat22 C = RotAT * RotB;
	Mat22 absC = Abs(C);
	Mat22 absCT = absC.Transpose();
	// Box A faces
	Vec2 faceA = Abs(dA) - hA - absC * hB;
	if (faceA.x > 0.0f || faceA.y > 0.0f)
		return 0;

	// Box B faces
	Vec2 faceB = Abs(dB) - absCT * hA - hB;
	if (faceB.x > 0.0f || faceB.y > 0.0f)
		return 0;

	// Find best axis
	Axis axis;
	float separation;
	Vec2 normal;

	// Box A faces
	axis = FACE_A_X;
	separation = faceA.x;
	normal = dA.x > 0.0f ? RotA.col1 : -RotA.col1;

	const float relativeTol = 0.95f;
	const float absoluteTol = 0.01f;

	if (faceA.y > relativeTol * separation + absoluteTol * hA.y)
	{
		axis = FACE_A_Y;
		separation = faceA.y;
		normal = dA.y > 0.0f ? RotA.col2 : -RotA.col2;
	}

	// Box B faces
	if (faceB.x > relativeTol * separation + absoluteTol * hB.x)
	{
		axis = FACE_B_X;
		separation = faceB.x;
		normal = dB.x > 0.0f ? RotB.col1 : -RotB.col1;
	}

	if (faceB.y > relativeTol * separation + absoluteTol * hB.y)
	{
		axis = FACE_B_Y;
		separation = faceB.y;
		normal = dB.y > 0.0f ? RotB.col2 : -RotB.col2;
	}

	// Setup clipping plane data based on the separating axis
	Vec2 frontNormal, sideNormal;
	ClipVertex incidentEdge[2];
	float front, negSide, posSide;
	char negEdge, posEdge;

	// Compute the clipping lines and the line segment to be clipped.
	switch (axis)
	{
	case FACE_A_X:
	{
		frontNormal = normal;
		front = Dot(posA, frontNormal) + hA.x;
		sideNormal = RotA.col2;
		float side = Dot(posA, sideNormal);
		negSide = -side + hA.y;
		posSide = side + hA.y;
		negEdge = EDGE3;
		posEdge = EDGE1;
		ComputeIncidentEdge(incidentEdge, hB, posB, RotB, frontNormal);
	}
	break;

	case FACE_A_Y:
	{
		frontNormal = normal;
		front = Dot(posA, frontNormal) + hA.y;
		sideNormal = RotA.col1;
		float side = Dot(posA, sideNormal);
		negSide = -side + hA.x;
		posSide = side + hA.x;
		negEdge = EDGE2;
		posEdge = EDGE4;
		ComputeIncidentEdge(incidentEdge, hB, posB, RotB, frontNormal);
	}
	break;

	case FACE_B_X:
	{
		frontNormal = -normal;
		front = Dot(posB, frontNormal) + hB.x;
		sideNormal = RotB.col2;
		float side = Dot(posB, sideNormal);
		negSide = -side + hB.y;
		posSide = side + hB.y;
		negEdge = EDGE3;
		posEdge = EDGE1;
		ComputeIncidentEdge(incidentEdge, hA, posA, RotA, frontNormal);
	}
	break;

	case FACE_B_Y:
	{
		frontNormal = -normal;
		front = Dot(posB, frontNormal) + hB.y;
		sideNormal = RotB.col1;
		float side = Dot(posB, sideNormal);
		negSide = -side + hB.x;
		posSide = side + hB.x;
		negEdge = EDGE2;
		posEdge = EDGE4;
		ComputeIncidentEdge(incidentEdge, hA, posA, RotA, frontNormal);
	}
	break;
	}

	// clip other face with 5 box planes (1 face plane, 4 edge planes)

	ClipVertex clipPoints1[2];
	ClipVertex clipPoints2[2];
	int np;

	// Clip to box side 1
	np = ClipSegmentToLine(clipPoints1, incidentEdge, -sideNormal, negSide, negEdge);

	if (np < 2)
		return 0;

	// Clip to negative box side 1
	np = ClipSegmentToLine(clipPoints2, clipPoints1, sideNormal, posSide, posEdge);

	if (np < 2)
		return 0;

	// Now clipPoints2 contains the clipping points.
	// Due to roundoff, it is possible that clipping removes all points.

	int numContacts = 0;
	for (int i = 0; i < 2; ++i)
	{
		float separation = Dot(frontNormal, clipPoints2[i].v) - front;

		if (separation <= 0)
		{
			contacts[numContacts].separation = separation;
			contacts[numContacts].normal = normal;
			// slide contact point onto reference face (easy to cull)
			contacts[numContacts].position = clipPoints2[i].v - separation * frontNormal;
			contacts[numContacts].feature = clipPoints2[i].fp;
			if (axis == FACE_B_X || axis == FACE_B_Y)
				Flip(contacts[numContacts].feature);
			++numContacts;
		}
	}

	return numContacts;
}

int CNCCollide(Contact*contacts, Body* bodyA, Body* bodyB)
{
	Vec2 posA = bodyA->position;
	Vec2 posB = bodyB->position;

	// dp = B물체의 중점에서 A물체의 중점으로 방향을 갖는 벡터
	Vec2 dp = posB - posA;

	//각 원들의 반지름
	float rA = bodyA->width.x * 0.5f;
	float rB = bodyB->width.x * 0.5f;

	//두 물체간에 거리
	float distSqr = Dot(dp, dp);
	float radius = rA + rB;

	//충돌 체크
	if (distSqr > radius * radius)
		return 0;

	Vec2 normal = posB - posA;
	normal.Normalize();
	
	Vec2 cA = posA + rA * normal;
	Vec2 cB = posB - rB * normal;
	int numContacts = 0;
	ClipVertex clipPoint;

	contacts[numContacts].separation = Dot(cB-cA, normal);
	contacts[numContacts].normal = normal;
	contacts[numContacts].position = 0.5f * (cA + cB);
	contacts[numContacts].feature = clipPoint.fp;
	++numContacts;
	bodyA->isItem = true;
	bodyB->isItem = true;
	return numContacts;
}

int CNBCollide(Contact*contacts, Body* bodyA, Body* bodyB)
{
	Vec2 posA = bodyA->position;
	Vec2 posB = bodyB->position;

	float rA = bodyA->width.x * 0.5f;
	// dp = B물체의 중점에서 A물체의 중점으로 방향을 갖는 벡터
	Vec2 dp = posB - posA;

	//네모의 모서리 길이 계산
	Vec2 d;
	d.Set(bodyB->width.x * 0.5f, bodyB->width.y * 0.5f);
	float Bx = bodyB->width.x * 0.5f;
	float By = bodyB->width.y * 0.5f;
	float rB = sqrtf(Bx * Bx + By * By);

	// Setup
	Vec2 hA = 0.5f * bodyA->width;
	Vec2 hB = 0.5f * bodyB->width;

	Mat22 RotA(bodyA->rotation), RotB(bodyB->rotation);

	Mat22 RotAT = RotA.Transpose();
	Mat22 RotBT = RotB.Transpose();

	Vec2 dB = RotBT * dp;
	Mat22 C = RotAT * RotB;
	Mat22 absC = Abs(C);
	Mat22 absCT = absC.Transpose();
	// Box A faces
	Vec2 faceB = Abs(dB) - hB - absC * hA;
	if (faceB.x > 0.0f || faceB.y > 0.0f)
		return 0;

	// Find best axis
	Axis axis;
	float separation;
	Vec2 normal;

	// Box A faces
	axis = FACE_B_X;
	separation = faceB.x;
	normal = dB.x > 0.0f ? RotB.col1 : -RotB.col1;

	const float relativeTol = 0.95f;
	const float absoluteTol = 0.01f;

	if (faceB.y > relativeTol * separation + absoluteTol * hB.y)
	{
		axis = FACE_B_Y;
		separation = faceB.y;
		normal = dB.y > 0.0f ? RotB.col2 : -RotB.col2;
	}

	int numContacts = 0;
	Vec2 frontNormal, sideNormal, cA, cB;
	ClipVertex incidentEdge[2];
	float front, negSide, posSide;
	char negEdge, posEdge;

	if (Dot(dp, dp) >= 1.0f)
	{
		//두 물체간에 거리
		float distSqr = Dot(dp, dp);
		float radius = rA + rB;

		//충돌 체크
		if (distSqr > radius * radius)
			return 0;

		Vec2 normal = posB - posA;
		normal.Normalize();

		Vec2 cA = posA + rA * normal;
		Vec2 cB = posB - rB * normal;

		contacts[numContacts].separation = Dot(cA - cB, normal);
		contacts[numContacts].normal = normal;
		contacts[numContacts].position = 0.5f * (cA + cB);
		++numContacts;
		return numContacts;
	}

	switch (axis)
	{
	case FACE_B_X:
	{
		frontNormal = normal;
		front = Dot(posB, frontNormal) + hB.x;
		sideNormal = RotB.col2;
		float side = Dot(posB, sideNormal);
		negSide = -side + hB.y;
		posSide = side + hB.y;
		negEdge = EDGE3;
		posEdge = EDGE1;
		ComputeIncidentEdge(incidentEdge, hA, posA, RotB, frontNormal);

		cB = posA + bodyA->width.x * 0.5f * frontNormal;
		cA = posB - rB * frontNormal;

	}
	break;

	case FACE_B_Y:
	{
		frontNormal = normal;
		front = Dot(posB, frontNormal) + hB.y;
		sideNormal = RotB.col1;
		float side = Dot(posB, sideNormal);
		negSide = -side + hB.x;
		posSide = side + hB.x;
		negEdge = EDGE2;
		posEdge = EDGE4;
		ComputeIncidentEdge(incidentEdge, hA, posA, RotB, frontNormal);

		cB = posA + bodyA->width.y * 0.5f * frontNormal;
		cA = posB - rB * frontNormal;
	}
	break;
	}
	ClipVertex clipPoints1[2];
	ClipVertex clipPoints2[2];
	int np;

	// Clip to box side 1
	np = ClipSegmentToLine(clipPoints1, incidentEdge, -sideNormal, negSide, negEdge);

	if (np < 2)
		return 0;

	// Clip to negative box side 1
	np = ClipSegmentToLine(clipPoints2, clipPoints1, sideNormal, posSide, posEdge);

	if (np < 2)
		return 0;

	contacts[numContacts].separation = Dot(cB - cA, frontNormal);
	contacts[numContacts].normal = frontNormal;
	contacts[numContacts].position = 0.5f * (clipPoints2[0].v - contacts[numContacts].separation * frontNormal + clipPoints2[1].v - contacts[numContacts].separation * frontNormal);
	contacts[numContacts].feature = clipPoints2[0].fp;
	++numContacts;
	bodyA->isJump = true;
	return numContacts;
}

int BNCCollide(Contact*contacts, Body* bodyA, Body* bodyB)
{
	Vec2 posA = bodyA->position;
	Vec2 posB = bodyB->position;

	float rB = bodyB->width.x * 0.5f;
	// dp = B물체의 중점에서 A물체의 중점으로 방향을 갖는 벡터
	Vec2 dp = posB - posA;

	//네모의 모서리 길이 계산
	Vec2 d;
	d.Set(bodyA->width.x * 0.5f, bodyA->width.y * 0.5f);
	float Ax = bodyA->width.x * 0.5f;
	float Ay = bodyA->width.y * 0.5f;
	float rA = sqrtf(Ax * Ax + Ay * Ay);

	// Setup
	Vec2 hA = 0.5f * bodyA->width;
	Vec2 hB = 0.5f * bodyB->width;

	Mat22 RotA(bodyA->rotation), RotB(bodyB->rotation);

	Mat22 RotAT = RotA.Transpose();
	Mat22 RotBT = RotB.Transpose();

	Vec2 dA = RotAT * dp;

	Mat22 C = RotAT * RotB;
	Mat22 absC = Abs(C);
	Mat22 absCT = absC.Transpose();
	// Box A faces
	Vec2 faceA = Abs(dA) - hA - absC * hB;
	if (faceA.x > 0.0f || faceA.y > 0.0f)
		return 0;

	// Find best axis
	Axis axis;
	float separation;
	Vec2 normal;

	// Box A faces
	axis = FACE_A_X;
	separation = faceA.x;
	normal = dA.x > 0.0f ? RotA.col1 : -RotA.col1;

	const float relativeTol = 0.95f;
	const float absoluteTol = 0.01f;

	if (faceA.y > relativeTol * separation + absoluteTol * hA.y)
	{
		axis = FACE_A_Y;
		separation = faceA.y;
		normal = dA.y > 0.0f ? RotA.col2 : -RotA.col2;
	}

	int numContacts = 0;
	Vec2 frontNormal, sideNormal, cA,cB;
	ClipVertex incidentEdge[2];
	float front, negSide, posSide;
	char negEdge, posEdge;

	if (Dot(dp, d) >= 1.0f)
	{
		//두 물체간에 거리
		float distSqr = Dot(dp, dp);
		float radius = rA + rB;
		//충돌 체크
		if (distSqr > radius * radius)
			return 0;

		Vec2 normal = posB - posA;
		normal.Normalize();

		Vec2 cA = posA + rA * normal;
		Vec2 cB = posB - rB * normal;

		contacts[numContacts].separation = Dot(cB - cA, normal);
		contacts[numContacts].normal = normal;
		contacts[numContacts].position = 0.5f * (cA + cB);
		++numContacts;
		return numContacts;
	}

	switch (axis)
	{
		case FACE_A_X:
		{
			frontNormal = normal;
			front = Dot(posA, frontNormal) + hA.x;
			sideNormal = RotA.col2;
			float side = Dot(posA, sideNormal);
			negSide = -side + hA.y;
			posSide = side + hA.y;
			negEdge = EDGE3;
			posEdge = EDGE1;
			ComputeIncidentEdge(incidentEdge, hB, posB, RotA, frontNormal);

			cA = posA + bodyA->width.x * 0.5f * frontNormal;
			cB = posB - rB * frontNormal;

		}
		break;

		case FACE_A_Y:
		{
			frontNormal = normal;
			front = Dot(posA, frontNormal) + hA.y;
			sideNormal = RotA.col1;
			float side = Dot(posA, sideNormal);
			negSide = -side + hA.x;
			posSide = side + hA.x;
			negEdge = EDGE2;
			posEdge = EDGE4;
			ComputeIncidentEdge(incidentEdge, hB, posB, RotA, frontNormal);

			cA = posA + bodyA->width.y * 0.5f * frontNormal;
			cB = posB - rB * frontNormal;
		}
		break;
	}
	ClipVertex clipPoints1[2];
	ClipVertex clipPoints2[2];
	int np;

	// Clip to box side 1
	np = ClipSegmentToLine(clipPoints1, incidentEdge, -sideNormal, negSide, negEdge);
	
	if (np < 2)
		return 0;

	// Clip to negative box side 1
	np = ClipSegmentToLine(clipPoints2, clipPoints1, sideNormal, posSide, posEdge);

	if (np < 2)
		return 0;

		contacts[numContacts].separation = Dot(cB - cA, frontNormal);
		contacts[numContacts].normal = frontNormal;
		contacts[numContacts].position = 0.5f * (clipPoints2[0].v - contacts[numContacts].separation * frontNormal + clipPoints2[1].v - contacts[numContacts].separation * frontNormal);
		contacts[numContacts].feature = clipPoints2[0].fp;
		++numContacts;
		bodyB->isJump = true;
	return numContacts;
}
// The normal points from A to B
int Collide(Contact* contacts, Body* bodyA, Body* bodyB)
{
	if (bodyA->shape == Body::box && bodyB->shape == Body::box)
		return BNBCollide(contacts, bodyA, bodyB);
	else if (bodyA->shape == Body::circle && bodyB->shape == Body::circle)
		return CNCCollide(contacts, bodyA, bodyB);
	else if (bodyA->shape == Body::box && bodyB->shape == Body::circle)
		return BNCCollide(contacts, bodyA, bodyB);
	else if (bodyA->shape == Body::circle && bodyB->shape == Body::box)
		return BNCCollide(contacts, bodyA, bodyB);
	else if (bodyA->shape == Body::triangle && bodyB->shape == Body::triangle)
		return BNCCollide(contacts, bodyA, bodyB);
	else
		return BNBCollide(contacts, bodyA, bodyB);
}