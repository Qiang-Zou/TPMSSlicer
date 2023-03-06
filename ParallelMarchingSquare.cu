#include "cuda.h"
#include "cuda_runtime.h"
#include<device_launch_parameters.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ParallelMarchingSquare.h"
#include "call_cuda.h"
#define THREADS_PER_BLOCK			512
#define BLOCKS_PER_GRID				256

extern __global__ void lookupTable(int* points, ContourEdge* edges, int* globalEid, int nodeNum, int3 imageRes, SCALAR cellSize, SCALAR* left, SCALAR height,int* Map);
extern __global__ void generateContour(ContourEdge* edges, int* globalEid, int* Map);



__device__
int lineTable[16][4] = {
		{-1, -1, -1, -1},
		{ 3,  0, -1, -1},
		{ 0,  1, -1, -1},
		{ 3,  1, -1, -1},
		{ 1,  2, -1, -1},
		{ 3,  0,  1,  2},
		{ 0,  2, -1, -1},
		{ 3,  2, -1, -1},
		{ 2,  3, -1, -1},
		{ 2,  0, -1, -1},
		{ 0,  1,  2,  3},
		{ 2,  1, -1, -1},
		{ 1,  3, -1, -1},
		{ 1,  0, -1, -1},
		{ 0,  3, -1, -1},
		{-1, -1, -1, -1},
};

__device__
int vertexTable[16][3] = {
	{-1,-1,-1},
    {1,2,3},
    {0,2,3},
    {2,3,-1},
    {0,1,3},
    {1,3,-1},
    {0,3,-1},
    {3,-1,-1},
    {0,1,2},
    {1,2,-1},
    {0,2,-1},
    {2,-1,-1},
    {0,1,-1},
    {1,-1,-1},
    {0,-1,-1},
    {-1,-1,-1},
};
extern "C" void call_lookupTable(int* points, ContourEdge* edges, int* globalEid, int nodeNum, int3 imageRes, SCALAR cellSize, SCALAR* left, SCALAR height,int* Map)
{
	lookupTable << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (points, edges, globalEid, nodeNum, imageRes,cellSize,left,height,Map);
}
__global__ void lookupTable(int* points, ContourEdge* edges,int* globalEid,int nodeNum, int3 imageRes, SCALAR cellSize, SCALAR* left, SCALAR height, int* Map)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int x, y;
	int resX = imageRes.x; int resY = imageRes.z;
	int EresX = resX - 1; 
	while (index < nodeNum)
	{
		x = index % imageRes.x;		y = index / imageRes.x;
		if (x<(resX-1)&&y<(resY-1))
		{
			int verts[4];
			verts[0] = points[index];
			verts[1] = points[index + 1];
			verts[2] = points[index + resX + 1];
			verts[3] = points[index + resX];

			int type = 0;
			for (int i = 0; i < 4; i++)
			{
				type += verts[i] * pow((double)2, (double)i);
			}
			int* config = lineTable[type];
			int s;
			int eid;
			for (int i = 0; i < 4; i++)
			{
				s = config[i];
				if (s != -1)
				{
					int stickid;
					switch (s) {
					case 0:
						stickid = (index / resX) * 2 * EresX + index % resX;
						break;
					case 1:
						stickid = ((index + 1) / resX) * 2 * EresX + (index + 1) % resX + EresX;
						break;
					case 2:
						stickid = ((index + resX) / resX) * 2 * EresX + (index + resX) % resX;
						break;
					case 3:
						stickid = (index / resX) * 2 * EresX + index % resX + EresX;
						break;

					}
					if (i == 0 || i == 2)
					{
						eid = atomicAdd(globalEid, 1);
						//printf("%d ", eid);
						edges[eid].id = eid;
						edges[eid].st = stickid;
						Map[stickid] = eid;
						//计算起始点坐标
						int horiStickId; bool isHoriStick = true;
						if ((edges[eid].st / EresX) % 2 == 0)
							horiStickId = edges[eid].st;
						else
						{
							horiStickId = edges[eid].st - EresX;
							isHoriStick = false;
						}
						int ptX = horiStickId % EresX; int ptY = horiStickId / EresX / 2;
						int ptX2,ptY2;
						if (isHoriStick)
						{
							ptX2 = ptX + 1; ptY2 = ptY;
						}
						else
						{
							ptX2 = ptX; ptY2 = ptY + 1;
						}
						edges[eid].stpt[0] = left[0] + cellSize * (ptX + ptX2) / 2.0;
						edges[eid].stpt[1] = left[1] + cellSize * (ptY + ptY2) / 2.0;
						edges[eid].stpt[2] = height;
						//着色顶点
						int* vconfig = vertexTable[type]; int v;
						for (int j = 0; j < 3; j++)
						{
							v = vconfig[j];
							if (v != -1)
							{
								switch (v) {
								case 0:
									edges[eid].v[j] = index;
									break;
								case 1:
									edges[eid].v[j] = index + 1;
									break;
								case 2:
									edges[eid].v[j] = index + resX + 1;
									break;
								case 3:
									edges[eid].v[j] = index + resX;
									break;
								}
							}
						}
					}
					else
					{
						edges[eid].ed = stickid;
						//计算终止点坐标
						int horiStickId; bool isHoriStick = true;
						if ((edges[eid].ed / EresX) % 2 == 0)
							horiStickId = edges[eid].ed;
						else
						{
							horiStickId = edges[eid].ed - EresX;
							isHoriStick = false;
						}
						int ptX = horiStickId % EresX; int ptY = horiStickId / EresX / 2;
						int ptX2, ptY2;
						if (isHoriStick)
						{
							ptX2 = ptX + 1; ptY2 = ptY;
						}
						else
						{
							ptX2 = ptX; ptY2 = ptY + 1;
						}
						edges[eid].edpt[0] = left[0] + cellSize * (ptX + ptX2) / 2.0;
						edges[eid].edpt[1] = left[1] + cellSize * (ptY + ptY2) / 2.0;
						edges[eid].edpt[2] = height;
						
					}
				}
			}
		}
		index += blockDim.x * gridDim.x;
	}

}

extern "C" void call_generateContour(ContourEdge* edges, int* globalEid,int* Map)
{
	generateContour << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (edges, globalEid,Map);
}
__global__ void generateContour(ContourEdge* edges, int* globalEid, int* Map)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	int edgeNum = *globalEid;

	while (index < edgeNum)
	{
		ContourEdge* e = &edges[index]; 
		int end = e->ed;
		e->next = Map[end];
		/*for (int i = 0; i < edgeNum; i++)
		{
			if (edges[i].st == e->ed)
			{
				e->next = edges[i].id;
				break;
			}
				
		}*/
		index += blockDim.x * gridDim.x;
	}
}