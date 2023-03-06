#pragma once
typedef struct edge
{
	int id, st, ed, next;
	float stpt[3] = { 0,0,0 }, edpt[3] = {0,0,0};
	int v[3] = { -1,-1,-1 };
	//int lpedgeNum=0;
} ContourEdge;