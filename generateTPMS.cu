#include "cuda.h"
#include "cutil.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "LDNIcudaSolid.h"
#include "call_cuda.h"
#include <Eigen/Dense>


extern __global__ void generateTPMS(bool* digitImage, double cellSize, int nodeNum, int3 imgRes, double left0, double left1, Eigen::Vector3d left, Eigen::Vector3d right, double layerLevel, int idx, int up, int low);

extern "C" void call_generateTPMS(bool* digitImage, double cellSize, int nodeNum, int3 imgRes, double left0, double left1, Eigen::Vector3d left, Eigen::Vector3d right, double layerLevel, int idx, int up, int low)
{
	generateTPMS << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (digitImage,cellSize,nodeNum,imgRes,left0,left1,left,right,layerLevel,idx,up,low);
}

__global__ void generateTPMS(bool* digitImage,double cellSize,int nodeNum, int3 imgRes, double left0,double left1,Eigen::Vector3d left, Eigen::Vector3d right,double layerLevel,int idx,int up,int low)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int rowIdx, colIdx;
	const double PI = 3.1415926;

	while (index < nodeNum)
	{
		rowIdx = index % imgRes.x;
		colIdx = index / imgRes.x;

		double x = left0 + rowIdx * cellSize;
		double y = left1 + colIdx * cellSize;
		double z = layerLevel;

		bool flag = !(x < left[0] || x > right[0] || y < left[1] || y > right[1] || z < left[2] || z > right[2]);
		if (!flag)
		{
			index += blockDim.x * gridDim.x;
			continue;
		}

		double lv;
		if (idx == 1)
		{
			lv = cos(x) + cos(y) + cos(z);
		}
		else if (idx == 2)
		{
			double l1(10), l2(3);
			double w1(2 * PI / l1), w2(2 * PI / l2);
			auto phiG = cos(w1 * x) * sin(w1 * y) + cos(w1 * y) * sin(w1 * z) + cos(w1 * z) * sin(w1 * x);
			auto phiGSon = cos(w2 * x) * sin(w2 * y) + cos(w2 * y) * sin(w2 * z) + cos(w2 * z) * sin(w2 * x);
			lv = min(phiG, phiGSon);
		}
		else if (idx == 3)
		{
			double l1(24), l2(6), l3(1.5);
			double w1(2 * PI / l1), w2(2 * PI / l2), w3(2 * PI / l3);
			auto Phi_G = cos(w1*x)*sin(w1*y) + cos(w1*y)*sin(w1*z) + cos(w1*z)*sin(w1*x) + 0.7;
			auto Phi_G2 = cos(w1*x)*sin(w1*y) + cos(w1*y)*sin(w1*z) + cos(w1*z)*sin(w1*x) - 0.7;
			auto Phi_G_son = cos(w2*x)*sin(w2*y) + cos(w2*y)*sin(w2*z) + cos(w2*z)*sin(w2*x) + 0.7;
			auto Phi_G_son_2 = cos(w2*x)*sin(w2*y) + cos(w2*y)*sin(w2*z) + cos(w2*z)*sin(w2*x) - 0.7;
			auto Phi_G_son_3 = cos(w3*x)*sin(w3*y) + cos(w3*y)*sin(w3*z) + cos(w3*z)*sin(w3*x) + 0.5;
			auto Phi_G_son_3_2 = cos(w3*x)*sin(w3*y) + cos(w3*y)*sin(w3*z) + cos(w3*z)*sin(w3*x) - 0.5;
			// first and second level;
			auto Phi1 = min(Phi_G, -Phi_G2);
			auto Phi2 = min(Phi_G_son, -Phi_G_son_2);
			auto Phi_f = min(Phi1, Phi2);
			// second and third level;
			auto Phi3 = min(Phi_G_son, -Phi_G_son_2);
			auto Phi4 = min(Phi_G_son_3, -Phi_G_son_3_2);
			auto Phi_son = min(Phi3, Phi4);
			lv = min(Phi_f, Phi_son);
		}
		else if (idx == 4)
		{
			double l1(24), l2(6), l3(1.5);
			double w1(2 * PI / l1), w2(2 * PI / l2), w3(2 * PI / l3);
			auto Phi_G = cos(w1*x)*sin(w1*y) + cos(w1*y)*sin(w1*z) + cos(w1*z)*sin(w1*x) + 0.7;
			auto Phi_G2 = cos(w1*x)*sin(w1*y) + cos(w1*y)*sin(w1*z) + cos(w1*z)*sin(w1*x) - 0.7;
			auto Phi_G_son = cos(w2*x)*sin(w2*y) + cos(w2*y)*sin(w2*z) + cos(w2*z)*sin(w2*x) + 0.7;
			auto Phi_G_son_2 = cos(w2*x)*sin(w2*y) + cos(w2*y)*sin(w2*z) + cos(w2*z)*sin(w2*x) - 0.7;
			auto Phi_G_son_3 = cos(w3*x)*sin(w3*y) + cos(w3*y)*sin(w3*z) + cos(w3*z)*sin(w3*x) + 0.5;
			auto Phi_G_son_3_2 = cos(w3*x)*sin(w3*y) + cos(w3*y)*sin(w3*z) + cos(w3*z)*sin(w3*x) - 0.5;
			// first and second level;
			auto Phi1 = min(Phi_G, -Phi_G2);
			auto Phi2 = min(Phi_G_son, -Phi_G_son_2);
			auto Phi_f = min(Phi1, Phi2);
			// second and third level;
			auto Phi3 = min(Phi_G_son, -Phi_G_son_2);
			auto Phi4 = min(Phi_G_son_3, -Phi_G_son_3_2);
			auto Phi_son = min(Phi3, Phi4);
			lv = min(Phi_f, Phi_son);
		}
		else
		{
			printf("Unknown type, not implemented yet!\n");
		}
		
		if (lv <= up && lv >= low)digitImage[rowIdx+ colIdx * imgRes.x] = true;

		index += blockDim.x * gridDim.x;
	}

}