#include "TPMSModeler.h"
#include "Utils/MarchingSquare/MarchingSquareBin.h"
#include "Utils/MarchingCube/MC33.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
#include <experimental/unordered_map>
#include <time.h>

#include <QImage>

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/remove.h>
#include <thrust/sequence.h>
#include <thrust/fill.h>
#include <thrust/count.h>
#include <thrust/functional.h>
#include <thrust/sequence.h>
#include <thrust/iterator/reverse_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/unique.h>
#include <cstdlib>
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cutil.h"
#include "LDNIcudaOperation.h"
#include "LDNIcudaSolid.h"
#include "call_cuda.h"
#define BLOCKSIZE	64

TPMSModeler::TPMSModeler() :
    tolerence(1e-6),
    pathOfFiles("Results/")
{
}

TPMSModeler::~TPMSModeler()
{
    // delete files
    // if(remove("../results/in_surface.off") != 0)
    //    std::cerr << "TPMSModeler > Error deleting files" << endl;
}

void TPMSModeler::readTPMSParameters()
{
    std::cout << "\nPlease select the surface type for slicing: 1. Schwarz-Primitive, 2. Multiscale, 3. Multiscale2, 4. TPMS+Dist. Field ";
    int idx; std::cin >> idx;

    double up, low;
    std::cout << "\nPlease specify upper and lower bounds for TPMS: ";
    std::cin >> up >> low;
    Eigen::Vector3d left, right;
    std::cout << "Please specify the bounding box of TPMS: ";
    std::cin >> left[0] >> left[1] >> left[2];
    std::cin >> right[0] >> right[1] >> right[2];
    //    int idx(4);
    //    double up(1000), low(0);
    //    Vector3d left(0, 0, 0), right(63, 42, 85);

    switch (idx) {
    case 1: ptrModel = std::static_pointer_cast<TPMS>(std::make_shared<TPMSSchwarz>(idx,up, low, left, right));
        break;
    case 2: ptrModel = std::static_pointer_cast<TPMS>(std::make_shared<TPMSMultiscale>(idx,up, low, left, right));
        break;
    case 3: ptrModel = std::static_pointer_cast<TPMS>(std::make_shared<TPMSMultiscale2>(idx,up, low, left, right));
        break;
    case 4: ptrModel = std::static_pointer_cast<TPMS>(std::make_shared<TPMSMultiscaleDistFieldBound>(idx,up, low, left, right, pathOfFiles + "bone.txt"));
        break;
    default: std::cout << "Unknown type, not implemented yet!\n";
    }
}

void TPMSModeler::readTPMSParameters(QMeshPatch *&ptrMesh)
{
    readTPMSParameters();
    long time = clock();
    ptrMesh = doMarchingCube();
    std::cout << "TPMSModeler > Marching Cube Time (second) for "<< float(clock() - time) / CLOCKS_PER_SEC << endl;
}

bool TPMSModeler::sliceModelImplicit()
{
    // input bounding box, layer thickness, resolution for slicing
    double layerThick;
    std::cout << "Please specify the layer thickness for slicing: ";
    std::cin >> layerThick;

    size_t resX, resY;
    std::cout << "Please choose the resolution of the discretization: ";
    std::cin >> resX >> resY;

	// Support structure
	bool *gridNodes;
	//vector<vector<int>> binaryNodes[40];
	int imageSize[3] = { 0,0,0 };

	int*** binaryNodes;
	binaryNodes = new int**[100];
	
    //    double layerThick = 10;
    //    size_t resX(800), resY(800);
    Eigen::Vector3d left, right;
    ptrModel->getBBox(left, right);

	double sizeX = right[0] - left[0], sizeY = right[1] - left[1];
	left += Eigen::Vector3d(-sizeX / 20, -sizeY / 20, 0);
	right += Eigen::Vector3d(sizeX / 20, sizeY / 20, 0);

	double cellSize = min((right[0] - left[0]) / min(resX, resY), (right[1] - left[1]) / min(resX, resY));
	resX = (right[0] - left[0]) / cellSize;
	resY = (right[1] - left[1]) / cellSize;
	imageSize[0] = resX; //cout << (right[0] - left[0]) << endl;
	imageSize[2] = resY; //cout << (right[1] - left[1]) << endl;
	int nodeNum = imageSize[0] * imageSize[2];
	int3 imgRes = make_int3(imageSize[0], imageSize[1], imageSize[2]);

	for (int i = 0; i < 100; i++)
	{
		binaryNodes[i] = new int*[resX];
		for (int j = 0; j < resX; j++)
		{
			binaryNodes[i][j] = new int[resY];
		}
	}

    long time = clock();

	stack<map<int, int>> stk1;
	stack<int* >stk2;
	stack<int>stk3;
	stack<ContourEdge*>stk4;
    // a functor for slicing at each iteration
    // the image stores in-pixel as 1 (only these pixels are stored) and out-pixel as 0
    auto generateBinImage = [&](double layerLevel, int layerNum, string path){
        // padding of the image to handle the openness of TPMS at boundaries
        //Eigen::Vector3d leftLocal(left), rightLocal(right);
        //double sizeX = right[0] - left[0], sizeY = right[1] - left[1];
        //leftLocal += Eigen::Vector3d(-sizeX / 20, -sizeY / 20, 0);
        //rightLocal += Eigen::Vector3d(sizeX / 20, sizeY / 20, 0);
        //leftLocal[2] = rightLocal[2] = layerLevel;

        //// determine the cell size
        //double cellSize(std::min((rightLocal[0] - leftLocal[0]) / std::min(resX, resY), (rightLocal[1] - leftLocal[1]) / std::min(resX, resY)));
        //size_t xRes(std::ceil((rightLocal[0] - leftLocal[0]) / cellSize));
        //size_t yRes(std::ceil((rightLocal[1] - leftLocal[1]) / cellSize));
        //if (xRes < yRes) yRes = std::max(yRes, std::max(resX, resY));
        //else xRes = std::max(xRes, std::max(resX, resY));
        std::cout << "TPMSModeler > digital image resoltuion at "<< layerNum << "-th layer :" << resX << "*" << resY << endl;

		long binaryImageTime = clock();

        //std::unordered_map<std::string, double> digitImage; // string: "rowIdx+space+colIdx"
		bool* digitImage;
		CUDA_SAFE_CALL(cudaMalloc((void **)&digitImage, nodeNum * sizeof(bool)));
		CUDA_SAFE_CALL(cudaMemset((void*)digitImage, false, nodeNum * sizeof(bool)));

        // compute each pixel's binary value
        /*for (int rowIdx(0); rowIdx <= resX; ++rowIdx) {
            for (int colIdx(0); colIdx <= resY; ++colIdx) {
                double x = left[0] + rowIdx * cellSize;
                double y = left[1] + colIdx * cellSize;
                double z = layerLevel;
                if (ptrModel->isInSolid(x, y, z)) {
                    auto key = std::to_string(rowIdx) + " " + std::to_string(colIdx);
                    digitImage.insert(std::make_pair(key, 1));
                }
            }
        }
        if(digitImage.empty()) return;*/
		int idx; double up, low;
		ptrModel->getIdx(idx);
		ptrModel->getBounds(up, low);
		Eigen::Vector3d modelLeft, modelRight;
		ptrModel->getBBox(modelLeft, modelRight);

		call_generateTPMS(digitImage,cellSize,nodeNum,imgRes, left[0],left[1], modelLeft, modelRight,layerLevel,idx,up,low);
		bool* hostDigit = new bool[nodeNum];
		CUDA_SAFE_CALL(cudaMemcpy((void*)hostDigit, (void*)digitImage,nodeNum * sizeof(bool), cudaMemcpyDeviceToHost));

		int* DImage = new int[nodeNum];
		
        // save the bindary image to file
        int width = resX, height = resY;
        QImage img(width, height, QImage::Format_RGB32);
        for (size_t i(0); i < width; ++i) {
            for (size_t j(0); j < height; ++j) {
                //auto iter = digitImage.find(std::to_string(i) + " " + std::to_string(j));
				//if (iter != digitImage.end())
				if(hostDigit[i + j * imgRes.x])
				{
					img.setPixel(i, j, qRgb(0, 0, 0));
					binaryNodes[layerNum][i][j] = 1;
					DImage[i + j * resX] = 1;
				}
				else
				{
					img.setPixel(i, j, qRgb(255, 255, 255));
					binaryNodes[layerNum][i][j] = 0;
					DImage[i + j * resX] = 0;
				}
            }
        }
        string path_ = pathOfFiles + "layer_" + std::to_string(layerNum) + ".jpg";
        if (!img.save(QString::fromStdString(path_), "JPEG", 50)) std::cerr << "TPMSModeler > error writing picture at layer: " << layerNum << endl;
        else std::cout << "TPMSModeler > successfully writing picture at layer: " << layerNum << endl;

		std::cout << "TPMSModeler > Binary Image Time (micro second) for " << layerNum << "-th layer:" << double(clock() - binaryImageTime) / CLOCKS_PER_SEC << endl;
        // do contouring with optimization
        //double leftCorner[2] = {0, 0};
        //MarchingSquareBin msb(digitImage, xRes, yRes, cellSize, layerLevel - left[2], leftCorner); // binImage will be moved, not copy, into MarchingSquareBin
        //msb.doContouringWithOptim();
        ////msb.doContouring();
        ////path = pathFiles + "layer_" + std::to_string(layerNum) + ".off";
        ////msb.writeContours(path);
        //auto result = msb.getContours(false);
        //.......................
        // do things with "result", e.g., save to "path"

        //........................
		long time = clock();

		

		int* cudaDImage;
		CUDA_SAFE_CALL(cudaMalloc((void **)&cudaDImage, nodeNum * sizeof(int)));
		CUDA_SAFE_CALL(cudaMemcpy((void*)cudaDImage, (void*)DImage, nodeNum * sizeof(int), cudaMemcpyHostToDevice));
		ContourEdge* edges;
		CUDA_SAFE_CALL(cudaMalloc((void **)&edges, 2 * (resX - 1)*(resY - 1) * sizeof(ContourEdge)));
		int* globalEid;
		CUDA_SAFE_CALL(cudaMalloc((void **)&globalEid, sizeof(int)));
		CUDA_SAFE_CALL(cudaMemset((void*)globalEid, 0, sizeof(int)));
		int *Map;
		CUDA_SAFE_CALL(cudaMalloc((void **)&Map, 2000000 * sizeof(int)));
		CUDA_SAFE_CALL(cudaMemset((void*)Map, 0, 2000000 * sizeof(int)));


		SCALAR* cudaLeft;
		SCALAR* tmpLeft = new SCALAR[3];
		tmpLeft[0] = left[0]; tmpLeft[1] = left[1]; tmpLeft[2] = left[2];

		CUDA_SAFE_CALL(cudaMalloc((void **)&cudaLeft, 3 * sizeof(SCALAR)));
		CUDA_SAFE_CALL(cudaMemcpy((void*)cudaLeft, (void*)tmpLeft, 3 * sizeof(SCALAR), cudaMemcpyHostToDevice));
		call_lookupTable(cudaDImage, edges, globalEid, nodeNum, imgRes, cellSize, cudaLeft, layerLevel, Map);
		call_generateContour(edges, globalEid, Map);

		//检查Map
		/*int* hostMap = new int[20000];
		CUDA_SAFE_CALL(cudaMemcpy((void*)hostMap, (void*)Map, 20000 * sizeof(int), cudaMemcpyDeviceToHost));*/

		ContourEdge* hostEdges = new ContourEdge[2 * (resX - 1)*(resY - 1)];
		CUDA_SAFE_CALL(cudaMemcpy((void*)hostEdges, (void*)edges, 2 * (resX - 1)*(resY - 1) * sizeof(ContourEdge), cudaMemcpyDeviceToHost));
		int* hostEid = new int[1];
		CUDA_SAFE_CALL(cudaMemcpy((void*)hostEid, (void*)globalEid, sizeof(int), cudaMemcpyDeviceToHost));

        std::cout << "TPMSModeler > successfully generating contours at layer: " << layerNum << endl;
		std::cout << "TPMSModeler > MarchingSquare Time (micro second) for " << layerNum << "-th layer:" << double(clock() - time) / CLOCKS_PER_SEC << endl;

		bool* visit = new bool[*hostEid];
		fill(visit, visit + *hostEid, false);

		int *lpedgeNum = new int[2000]();
		int  loopNum = 0, pixelNum = *hostEid;
		map<int, int>mp;//保存每个loop对应的起始斜边的id
		QImage contourImg(width, height, QImage::Format_RGB32);
		for (int i = 0; i < *hostEid; i++)
		{
			if (!visit[i])
			{
				loopNum++;
				mp[loopNum] = i;
				lpedgeNum[loopNum]++;
				visit[i] = true;
				ContourEdge te = hostEdges[i];
				int firstId = te.id;
				int* config = te.v;
				for (int j = 0; j < 3; j++)
				{
					int v = config[j];
					if (v != -1)
					{
						int row = v % imgRes.x;
						int col = v / imgRes.x;
						contourImg.setPixel(row, col, qRgb(0, 0, 0));
					}
				}
				ContourEdge nextEdge = hostEdges[te.next];
				int nextId = nextEdge.id;
				while (nextId != firstId)
				{
					visit[nextId] = true; lpedgeNum[loopNum]++;
					int* nconfig = nextEdge.v;
					for (int j = 0; j < 3; j++)
					{
						int v = nconfig[j];
						if (v != -1)
						{
							int row = v % imgRes.x;
							int col = v / imgRes.x;
							contourImg.setPixel(row, col, qRgb(0, 0, 0));
						}
					}
					nextEdge = hostEdges[nextEdge.next];
					nextId = nextEdge.id;
				}
			}

		}
		cout << "loopNum: " << loopNum << endl;
		string offpath = pathOfFiles + "layer_" + std::to_string(layerNum) + ".contours.off";
		outputOFF(offpath, pixelNum, lpedgeNum, loopNum, mp, hostEdges);

		stk1.push(mp);
		stk2.push(lpedgeNum);
		stk3.push(*hostEid);
		stk4.push(hostEdges);

		std::string cpath = pathOfFiles + "layer_" + std::to_string(layerNum) + ".contour.jpg";
		if (!contourImg.save(QString::fromStdString(cpath), "JPEG", 50)) std::cerr << "Contour > error writing picture at layer: " << layerNum << endl;
		else std::cout << "Contour > successfully writing picture at layer: " << layerNum << endl;

    };

    // slice the model layer-by-layer
    double zMin(left[2]), zMax(right[2]);
    size_t layerNum(1);
    double layerLevel, layerBase;
    layerLevel = layerBase = zMin + layerThick;
    zMax -= layerThick / 10;


    // slicing loop
    string path = pathOfFiles + "contours" + ".cli";
    while (layerLevel < zMax) {
        // do slicing
        generateBinImage(layerLevel, layerNum, path);
        std::cout << "TPMSModeler > Slicing Time (micro second) for "<< layerNum << "-th layer:" << float(clock() - time) / CLOCKS_PER_SEC << endl;
        time = clock();

        ++layerNum;
        layerLevel += layerThick;
    }
    std::cout << "TPMSModeler > successfully slicing the model" << endl;

    //    //    // write to .cli files
    //    //    path = pathOfFiles + "contours1" + ".cli";
    //    //    writeCLIFileBin(path, layers);
    //    //    path = pathOfFiles + "contours1" + ".cli";
    //    //    writeCLIFile(path, layers);
    //    path = pathOfFiles + "contours" + ".cli";
    //    convertCLIFileBin2ASCII(path);

    //    std::cout << "TPMSModeler > successfully saving the contours" << endl;
	// write to .cli files
	int stkSize = stk1.size();
	for (int i = 0; i < stkSize; i++)
	{
		map<int, int> mp = stk1.top();
		mpvector.push_back(mp);
		stk1.pop();
		int* lpedgeNum = stk2.top();
		lpenumvector.push_back(lpedgeNum);
		stk2.pop();
		int eNum = stk3.top();
		edgeNum.push_back(eNum);
		stk3.pop();
		ContourEdge* edges = stk4.top();
		layer.push_back(edges);
		stk4.pop();

	}
	string clipath = "Results/CLIFileforPart/out_volume.contours1.cli";
	writeCLI(clipath, layer, edgeNum, mpvector, lpenumvector);

	long suptTime = clock();
	// Support structure
	imageSize[1] = layerNum - 1;
	CUDA_SAFE_CALL(cudaMalloc((void**)&(gridNodes), imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)gridNodes, 0.0, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool)));

	int idx = 0;
	bool* tmp = new bool[imageSize[0] * imageSize[1] * imageSize[2]];
	for (int j = 0; j < imageSize[2]; j++)
	{
		for (int k = 1; k <= imageSize[1]; k++)
		{

			for (int i = 0; i < imageSize[0]; i++)
			{
				//gridNodes[idx] = binaryNodes[i][j];
				if (binaryNodes[k][i][j] == 1)
					tmp[idx] = true;
				else
					tmp[idx] = false;
				idx++;
			}

		}
	}
	CUDA_SAFE_CALL(cudaMemcpy((void*)gridNodes, (void*)tmp, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool), cudaMemcpyHostToDevice));
	//User specified variables
	double anchorR = 2;
	double thres = 0.2;
	double nSampleWidth = 0.005;
	double cylinderRadius = 4.5;
	double patternThickness = 3;

	bool *suptNodes;
	bool *tempImage;
	bool *targetImage;
	bool *assistImage;
	bool *temp3D;


	double anchorRadius = anchorR;
	double threshold = thres;

	int suptRadius = (int)floor(anchorRadius*1.414 / nSampleWidth);
	int suptGridResX = (imageSize[0] - 1) / suptRadius + 1;
	int suptGridResZ = (imageSize[2] - 1) / suptRadius + 1;

	imgRes = make_int3(imageSize[0], imageSize[1], imageSize[2]);
	int3 suptimgRes = make_int3(imageSize[0], imageSize[1] - 1, imageSize[2]);
	//int nodeNum = imageSize[0] * imageSize[2];

	long ti[10] = { 0 };

	CUDA_SAFE_CALL(cudaMalloc((void**)&(temp3D), imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)temp3D, 0, imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));


	CUDA_SAFE_CALL(cudaMalloc((void**)&(suptNodes), imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)suptNodes, false, imageSize[0] * (imageSize[1] - 1)*imageSize[2] * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(tempImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)tempImage, false, nodeNum * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(targetImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)targetImage, false, nodeNum * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(assistImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)assistImage, false, nodeNum * sizeof(bool)));


	short2 **disTextures = (short2 **)malloc(2 * sizeof(short2 *));

	int disTexSize = max(imageSize[0], imageSize[2]);
	int factor = ceil((float)disTexSize / BLOCKSIZE);
	disTexSize = BLOCKSIZE * factor;

	int disMemSize = disTexSize * disTexSize * sizeof(short2);

	// Allocate 2 textures

	cudaMalloc((void **)&disTextures[0], disMemSize);
	cudaMalloc((void **)&disTextures[1], disMemSize);

	unsigned int *LinkIndex;
	CUDA_SAFE_CALL(cudaMalloc((void **)&LinkIndex, (nodeNum + 1) * sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemset((void*)LinkIndex, 0, (nodeNum + 1) * sizeof(unsigned int)));

	long t;

	for (int i = imageSize[1] - 2; i > -1; i--)
	{

		t = clock();
		call_krSLAContouring_Initialization(tempImage, targetImage, gridNodes, nodeNum, imgRes, i);
		CUDA_SAFE_CALL(cudaMemcpy(assistImage, targetImage, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));
		ti[0] += clock() - t;

		t = clock();
		LDNIcudaOperation::LDNISLAContouring_GrowthAndSwallow(threshold, targetImage, tempImage, i, imageSize, nSampleWidth, disTextures[0], disTextures[1], disTexSize);


		ti[1] += clock() - t;
		//add new support cylinder if necessary
		CUDA_SAFE_CALL(cudaMemset((void*)tempImage, false, nodeNum * sizeof(bool)));
		CUDA_SAFE_CALL(cudaMemset((void*)assistImage, false, nodeNum * sizeof(bool)));
		call_krSLAContouring_Filter1(assistImage, tempImage, targetImage, suptGridResX*suptGridResZ, make_int2(imageSize[0], imageSize[2]), suptRadius, i);


		//first step: support region growth in first class cylinders
		LDNIcudaOperation::LDNISLAContouring_GrowthAndSwallow(anchorRadius, targetImage, assistImage, i, imageSize, nSampleWidth, disTextures[0], disTextures[1], disTexSize);


		//second step: prepare second class cylinders and perform support region growth in second class cylinders
		CUDA_SAFE_CALL(cudaMemset((void*)assistImage, false, nodeNum * sizeof(bool)));
		call_krSLAContouring_OrthoSearchRemainAnchorZ(assistImage, tempImage, targetImage, suptGridResX,
			make_int2(suptRadius, imageSize[2]), make_int2(imageSize[0], imageSize[2]), i);



		call_krSLAContouring_OrthoSearchRemainAnchorX(assistImage, tempImage, targetImage, suptGridResZ,
			make_int2(suptRadius, imageSize[0]), make_int2(imageSize[0], imageSize[2]), i);



		LDNIcudaOperation::LDNISLAContouring_GrowthAndSwallow(anchorRadius, targetImage, assistImage, i, imageSize, nSampleWidth, disTextures[0], disTextures[1], disTexSize);
		long time = clock();

		//third step: prepare third class cylinders and support region growth in all third class cylinders
		CUDA_SAFE_CALL(cudaMemset((void*)assistImage, false, nodeNum * sizeof(bool)));
		LDNIcudaOperation::LDNISLAContouring_ThirdClassCylinder(anchorRadius, targetImage, assistImage, tempImage, make_int2(imageSize[0], imageSize[2]), nSampleWidth, i, disTextures[0], disTextures[1], disTexSize);


		//generate support structure map for this layer and update the cylinder position information arrays

		call_krSLAContouring_Filter5(gridNodes, tempImage, suptNodes, LinkIndex, nodeNum, imgRes, i);




		call_krFDMContouring_CopyNodesrom2Dto3D(tempImage, temp3D, nodeNum, suptimgRes, i);

		//std::cout << "Generate support map time for " << i + 1 << "-th layer:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
	}

	cudaFree(disTextures[0]);
	cudaFree(disTextures[1]);
	free(disTextures);
	cudaFree(assistImage);
	cudaFree(targetImage);
	cudaFree(tempImage);

	unsigned int LinkNum;
	call_func::setdev_ptr(LinkIndex, nodeNum, LinkNum);


	short *linkLayerD;
	short2 *linkID;
	unsigned int *linkLayerC;
	unsigned int *temp2D;

	CUDA_SAFE_CALL(cudaMalloc((void **)&temp2D, nodeNum * sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemset((void*)temp2D, 0, nodeNum * sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&linkLayerC, LinkNum * sizeof(unsigned int)));
	//CUDA_SAFE_CALL(cudaMemset((void*)linkLayerC, imgRes.y+1, LinkNum*sizeof(unsigned int) ) );
	CUDA_SAFE_CALL(cudaMalloc((void **)&linkLayerD, LinkNum * sizeof(short)));
	CUDA_SAFE_CALL(cudaMemset((void*)linkLayerD, 0, LinkNum * sizeof(short)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&linkID, LinkNum * sizeof(short2)));
	CUDA_SAFE_CALL(cudaMemset((void*)linkID, 0, LinkNum * sizeof(short2)));


	call_func::setdev_ptr2(linkLayerC, LinkNum, imgRes);

	time = clock();
	call_krSLAContouring_FindAllLinks(LinkIndex, linkLayerC, linkLayerD, linkID, temp3D, temp2D, suptimgRes.x*suptimgRes.y*suptimgRes.z, suptimgRes);

	cudaFree(temp3D);
	cudaFree(temp2D);



	call_krSLAContouring_RelateAllLinksBetweenLayers(LinkIndex, linkLayerC, gridNodes, suptNodes, suptimgRes.x*suptimgRes.y*suptimgRes.z, imgRes);

	cudaFree(LinkIndex);
	std::cout << "Relate all links time: " << double(clock() - time) / CLOCKS_PER_SEC << endl;



	CUDA_SAFE_CALL(cudaMalloc((void**)&(tempImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)tempImage, false, nodeNum * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(targetImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)targetImage, false, nodeNum * sizeof(bool)));


	double realThreshold = (cylinderRadius*nSampleWidth - nSampleWidth) / nSampleWidth;
	int gridRadius = (int)floor(realThreshold);
	time = clock();
	for (int i = imageSize[1] - 2; i > -1; i--)
	{

		call_krFDMContouring_CopyNodesrom3Dto2D(targetImage, suptNodes, nodeNum, suptimgRes, i);
		CUDA_SAFE_CALL(cudaMemcpy(tempImage, targetImage, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));
		call_krFDMContouring_Dilation(targetImage, tempImage, nodeNum, imgRes, realThreshold, gridRadius, i);


		call_krFDMContouring_CopyNodesrom2Dto3D(tempImage, suptNodes, nodeNum, suptimgRes, i);

	}
	std::cout << "Dilation time:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
	cudaFree(targetImage);
	cudaFree(tempImage);
	time = clock();
	for (int i = imageSize[1] - 3; i > -1; i--)
		call_krFDMContouring_VerticalSpptPxlProp(gridNodes, suptNodes, nodeNum, imgRes, i);
	std::cout << "VerticalSpptPxlProp time:" << double(clock() - time) / CLOCKS_PER_SEC << endl;

	int linkThreshold = (int)(0.4 / nSampleWidth);
	int lengthofLayer = linkThreshold / 8;
	int furtherStepLength = lengthofLayer / 2;



	time = clock();
	LDNIcudaOperation::LDNISLAContouring_GenerateConnectionforCylinders(linkLayerC, linkLayerD, linkID, gridNodes, suptNodes, imageSize, linkThreshold,
		lengthofLayer, furtherStepLength, LinkNum, nSampleWidth);
	std::cout << "Generate Connection for Cylinders time:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
	cudaFree(linkLayerC);
	cudaFree(linkLayerD);
	cudaFree(linkID);
	//cudaFree(gridNodes);


	realThreshold = (patternThickness*nSampleWidth - nSampleWidth) / nSampleWidth;
	gridRadius = (int)floor(realThreshold);

	CUDA_SAFE_CALL(cudaMalloc((void**)&(tempImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)tempImage, false, nodeNum * sizeof(bool)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(targetImage), nodeNum * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset((void*)targetImage, false, nodeNum * sizeof(bool)));
	time = clock();
	for (int i = imageSize[1] - 2; i > -1; i--)
	{

		call_krFDMContouring_CopyNodesrom3Dto2D(targetImage, suptNodes, nodeNum, suptimgRes, i);
		CUDA_SAFE_CALL(cudaMemcpy(tempImage, targetImage, nodeNum * sizeof(bool), cudaMemcpyDeviceToDevice));
		call_krFDMContouring_Dilation(targetImage, tempImage, nodeNum, imgRes, realThreshold, gridRadius, i);
		call_krFDMContouring_CopyNodesrom2Dto3D(tempImage, suptNodes, nodeNum, suptimgRes, i);

	}
	std::cout << "Dilation time:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
	cudaFree(targetImage);
	cudaFree(tempImage);


	bool* gridtmp = new bool[imageSize[0] * imageSize[1] * imageSize[2]];
	CUDA_SAFE_CALL(cudaMemcpy((void*)gridtmp, (void*)gridNodes, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(bool), cudaMemcpyDeviceToHost));
	int cnt = 0;
	size_t width = imageSize[0], height = imageSize[2], slicenum = imageSize[1];
	QImage *img = new QImage[slicenum];
	for (int i = 0; i < slicenum; i++)
	{
		QImage tmpimg(width, height, QImage::Format_RGB32);
		img[i] = tmpimg;
	}
	for (int j = 0; j < imageSize[2]; j++)
	{

		for (int k = 0; k < imageSize[1]; k++)
		{

			for (int i = 0; i < imageSize[0]; i++)
			{
				if (gridtmp[cnt++])img[k].setPixel(i, j, qRgb(0, 0, 0));
				else img[k].setPixel(i, j, qRgb(255, 255, 255));

			}
		}

	}
	for (int i = 1; i <= imageSize[1]; i++)
	{
		std::string path = "Results/Support/out_volume." + std::to_string(i) + ".origin.jpg";
		if (!img[i - 1].save(QString::fromStdString(path), "JPEG", 50)) std::cerr << "Origin > error writing picture at layer: " << i << endl;
		else std::cout << "Origin > successfully writing picture at layer: " << i << endl;
	}


	std::unordered_map<std::string, double> *digitImage = new std::unordered_map<std::string, double>[slicenum - 1];//generate contours for support structure
	std::vector<std::shared_ptr<QMeshPatch>> suptlayers;
	layerLevel = zMin + layerThick;
	bool* supttmp = new bool[imageSize[0] * (imageSize[1] - 1)*imageSize[2]];
	CUDA_SAFE_CALL(cudaMemcpy((void*)supttmp, (void*)suptNodes, imageSize[0] * (imageSize[1] - 1) * imageSize[2] * sizeof(bool), cudaMemcpyDeviceToHost));
	cnt = 0;
	QImage *suptimg = new QImage[slicenum - 1];
	for (int i = 0; i < slicenum - 1; i++)
	{
		QImage tmpimg(width, height, QImage::Format_RGB32);
		suptimg[i] = tmpimg;
	}
	for (int j = 0; j < imageSize[2]; j++)
	{
		for (int k = 0; k < imageSize[1] - 1; k++)
		{

			for (int i = 0; i < imageSize[0]; i++)
			{
				if (supttmp[cnt++])
				{
					suptimg[k].setPixel(i, j, qRgb(255, 0, 0));
					auto key = std::to_string(i) + " " + std::to_string(j);
					digitImage[k].insert(std::make_pair(key, 1));
				}
				else suptimg[k].setPixel(i, j, qRgb(255, 255, 255));

			}
		}

	}
	for (int i = 1; i <= imageSize[1] - 1; i++)
	{
		std::string path = "Results/Support/out_volume." + std::to_string(i) + ".support.jpg";
		if (!suptimg[i - 1].save(QString::fromStdString(path), "JPEG", 50)) std::cerr << "Support > error writing picture at layer: " << i << endl;
		else std::cout << "Support > successfully writing picture at layer: " << i << endl;

		//long time = clock();
		//MarchingSquareBin msb(digitImage[i - 1], resX, resY, cellSize, layerLevel, left.data());
		//msb.doContouring();
		//path = "Results/Support/out_volume." + std::to_string(i) + ".contours.off";
		////msb.writeContours(path);
		//auto result = msb.getContours(true);
		//suptlayers.push_back(std::shared_ptr<QMeshPatch>(result));
		//layerLevel += layerThick;
		//std::cout << "Support >  MarchingSquare Time (micro second) for " << i << "-th layer:" << double(clock() - time) / CLOCKS_PER_SEC << endl;
	}
	// write to .cli files
	/*string suptclipath = "Results/Support/CLIFileforSupt/out_volume.contours.cli";
	writeCLIFileBin(suptclipath, suptlayers);
	suptclipath = "Results/Support/CLIFileforSupt/out_volume.contours1.cli";
	writeCLIFile(suptclipath, suptlayers);*/
	std::cout << "TPMSModeler > Support Generation Time (micro second):" << float(clock() - suptTime) / CLOCKS_PER_SEC << endl;
    return true;
}

QMeshPatch *TPMSModeler::doMarchingCube()
{
    // determine the cell size
    Eigen::Vector3d left, right;
    ptrModel->getBBox(left, right);

    double cellSize = 0.1; // change to a smaller value if you want a finer mesh
    size_t resX = std::ceil((right[0] - left[0]) / cellSize);
    size_t resY = std::ceil((right[1] - left[1]) / cellSize);
    size_t resZ = std::ceil((right[2] - left[2]) / cellSize);

    // populate the grid
    MCGrid g;
    g.set_grid_dimensions(resX + 1, resY + 1, resZ + 1);
    g.set_r0(left[0], left[1], left[2]); // the left corner
    g.set_ratio_aspect(cellSize, cellSize, cellSize);
    g.set_Ang(90, 90, 90);
    for (size_t i(0); i <= resX; ++i) {
        for (size_t j(0); j <= resY; ++j) {
            for (size_t k(0); k <= resZ; ++k) {
                double x(i * cellSize + left[0]), y(j * cellSize + left[1]), z(k * cellSize + left[2]);
                float value = ptrModel->getValue(x, y, z);
                g.set_grid_value(i, j, k, value);
            }
        }
    }

    // do the contouring
    float isovalue = 0;
    MCGen MC; MC.set_grid3d(&g);
    std::shared_ptr<MCSurface> s(MC.calculate_isosurface(isovalue));

    return nullptr;

    // convert to QMeshPatch and display
    QMeshPatch* mesh = new QMeshPatch;
    auto& nodeList = mesh->GetNodeList();
    auto& edgeList = mesh->GetEdgeList();
    auto& faceList = mesh->GetFaceList();
    GLKPOSITION Pos;
    GLKPOSITION PosNode;
    QMeshNode *node,*startNode,*endNode;
    QMeshEdge *edge;
    QMeshFace *face;
    vector<QMeshNode *> nodeArray;
    int faceNum,nodeNum,i, j;

    // handle vertices
    nodeNum = s->get_num_vertices();
    for(i=0;i<nodeNum;i++) {
        auto p = s->getVertex(i);
        node=new QMeshNode;
        node->SetMeshPatchPtr(mesh);
        node->SetCoord3D(p[0],p[1],p[2]);
        node->SetIndexNo(nodeList.GetCount()+1);
        nodeList.AddTail(node);
    }

    nodeArray.resize(nodeNum);
    i=0;
    for(Pos=nodeList.GetHeadPosition();Pos!=NULL;i++) {
        node=(QMeshNode*)(nodeList.GetNext(Pos));
        nodeArray[i]=node;
    }

    // handle faces
    faceNum = s->get_num_triangles();
    for(i=0;i<faceNum;i++) {
        int nodeIndex;
        face=new QMeshFace;
        face->SetMeshPatchPtr(mesh);
        face->SetIndexNo(faceList.GetCount()+1);
        faceList.AddTail(face);

        auto vertIdx = s->getTriangle(i);
        for (j = 0; j < 3; ++j) {
            nodeIndex = vertIdx[j];
            (face->GetAttachedList()).AddTail(nodeArray[nodeIndex]);
        }
    }

    // handle edges
    for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
        face=(QMeshFace*)(faceList.GetNext(Pos));

        int edgeNum=(face->GetAttachedList()).GetCount();
        face->SetEdgeNum(edgeNum);

        nodeArray.resize(edgeNum);
        i=0;
        for(PosNode=(face->GetAttachedList()).GetHeadPosition();PosNode!=NULL;i++) {
            nodeArray[i]=(QMeshNode*)((face->GetAttachedList()).GetNext(PosNode));
            (nodeArray[i]->GetFaceList()).AddTail(face);
        }
        for(i=0;i<edgeNum;i++) {
            edge=NULL;
            startNode=nodeArray[i];
            endNode=nodeArray[(i+1)%edgeNum];
            bool bDir;
            for(PosNode=(startNode->GetEdgeList()).GetHeadPosition();PosNode!=NULL;) {
                QMeshEdge *temp=(QMeshEdge *)((startNode->GetEdgeList()).GetNext(PosNode));
                if ((temp->GetStartPoint()==startNode) && (temp->GetEndPoint()==endNode)) {
                    edge=temp;	bDir=true;
                }
                else if ((temp->GetStartPoint()==endNode) && (temp->GetEndPoint()==startNode)) {
                    edge=temp;	bDir=false;
                }
            }
            if (edge) {
                face->SetEdgeRecordPtr(i,edge);
                if (bDir) {
                    face->SetDirectionFlag(i,true);
                    edge->SetLeftFace(face);
                }
                else {
                    face->SetDirectionFlag(i,false);
                    edge->SetRightFace(face);
                }
            }
            else {
                edge=new QMeshEdge;
                edge->SetMeshPatchPtr(mesh);
                edge->SetStartPoint(startNode);
                edge->SetEndPoint(endNode);
                edge->SetIndexNo(edgeList.GetCount()+1);
                edgeList.AddTail(edge);

                edge->SetLeftFace(face);
                face->SetEdgeRecordPtr(i,edge);
                face->SetDirectionFlag(i,true);
                (startNode->GetEdgeList()).AddTail(edge);
                (endNode->GetEdgeList()).AddTail(edge);
            }
        }

        face->GetAttachedList().RemoveAll();
    }

    //compute the normal
    for(Pos=faceList.GetHeadPosition();Pos!=NULL;) {
        face=(QMeshFace*)(faceList.GetNext(Pos));
        face->CalPlaneEquation();
    }
    return mesh;
}

void TPMSModeler::outputOFF(string &path, int pixelNum, int* lpedgeNum, int loopNum, map<int, int> &mp, ContourEdge* edges)
{
	char* filename = new char[strlen(path.c_str()) + 1];
	strcpy(filename, path.c_str());

	SCALAR x, y, z;

	// do writing
	std::ofstream outfile(filename);
	if (!outfile.is_open()) {
		printf("===============================================\n");
		printf("Can not open the data file - OFF File Export!\n");
		printf("===============================================\n");
		return;
	}
	// 1-2 line: node numbers edge numbers
	outfile << "OFF" << endl;
	outfile << pixelNum << " " << loopNum << " " << 0;
	// vertices
	for (int i = 0; i < pixelNum; i++)
	{
		x = edges[i].stpt[0];
		y = edges[i].stpt[1];
		z = edges[i].stpt[2];
		outfile << endl << x << " " << y << " " << z;
	}

	// faces (indices starting at 0)
	/*int num;
	auto t = faceList;
	for (auto iter = FaceBegin(); iter != FaceEnd(); iter = iter->Next()) {
		auto face = (QMeshFace*)(iter->data);
		num = face->GetEdgeNum();
		outfile << endl << num;
		for (i = 0; i < num; i++) outfile << " " << face->GetNodeRecordPtr(i)->GetIndexNo();
	}*/
	int num;
	for (int i = 1; i <= loopNum; i++)
	{
		num = lpedgeNum[i];
		outfile << endl << num;

		outfile << " " << edges[mp[i]].id;
		ContourEdge e = edges[mp[i]];
		int firstId = e.id;
		int nextId = e.next;
		while (nextId != firstId)
		{
			ContourEdge ne = edges[nextId];
			outfile << " " << ne.id;
			nextId = ne.next;
		}

	}
	outfile.close();
}

bool TPMSModeler::writeCLI(std::string &path, std::vector<ContourEdge*> &layer, vector<int> &edgeNum, vector<map<int, int>> &mpvector, vector<int*> &lpenumvector)
{
	// bounding box
	Eigen::Vector3d left(1e8 * Eigen::Vector3d::Ones()), right(-1e8 * Eigen::Vector3d::Ones()), p;
	int cnt = 0;
	for (ContourEdge*& l : layer) {
		int num = edgeNum[cnt++];
		for (int i = 0; i < num; i++) {
			p[0] = l[i].stpt[0];
			p[1] = l[i].stpt[1];
			p[2] = l[i].stpt[2];
			left[0] = std::min(left[0], p[0]);
			left[1] = std::min(left[1], p[1]);
			left[2] = std::min(left[2], p[2]);
			right[0] = std::max(right[0], p[0]);
			right[1] = std::max(right[1], p[1]);
			right[2] = std::max(right[2], p[2]);
		}
	}

	// shift the bottom to the xoy plane
	vector<double> layerHeights;
	Eigen::Vector3d offset(left[0] - 0.01, left[1] - 0.01, left[2] - 0.01);
	cnt = 0;
	for (ContourEdge*& l : layer) {
		int num = edgeNum[cnt++];
		for (int i = 0; i < num; i++) {
			p[0] = l[i].stpt[0];
			p[1] = l[i].stpt[1];
			p[2] = l[i].stpt[2];
			p[0] -= offset[0]; p[1] -= offset[1]; p[2] -= offset[2];
			l[i].stpt[0] = p[0];
			l[i].stpt[1] = p[1];
			l[i].stpt[2] = p[2];
		}

		layerHeights.push_back(p[2]);
	}
	left -= offset; right -= offset;

	// write to the file
	// a functor for determining ccw
	auto areaCal = [](ContourEdge* edges, int firstId) {
		double area(0);
		double p1[3], p2[3];
		ContourEdge firstEdge = edges[firstId];
		int nextId = firstEdge.next;
		while (nextId != firstId) {
			ContourEdge ne = edges[nextId];
			p1[0] = firstEdge.stpt[0]; p1[1] = firstEdge.stpt[1]; p1[2] = firstEdge.stpt[2];
			p2[0] = ne.stpt[0]; p2[1] = ne.stpt[1]; p2[2] = ne.stpt[2];
			area += p1[0] * p2[1] - p1[1] * p2[0];
			nextId = ne.next;
			firstEdge = ne;
		}
		ContourEdge ne = edges[firstId];
		p1[0] = firstEdge.stpt[0]; p1[1] = firstEdge.stpt[1]; p1[2] = firstEdge.stpt[2];
		p2[0] = ne.stpt[0]; p2[1] = ne.stpt[1]; p2[2] = ne.stpt[2];
		area += p1[0] * p2[1] - p1[1] * p2[0];

		return area / 2;
	};

	std::ofstream outfile;
	outfile.open(path, ios::out | ios::binary);
	if (!outfile.is_open()) {
		std::cerr << "LatticeModeler > error openning CLI file" << endl;
		return false;
	}
	// write header
	outfile << "$$HEADERSTART" << std::endl;
	outfile << "$$ASCII" << std::endl << "$$UNITS\/00000000.010000" << std::endl << "$$VERSION\/100" << std::endl;
	outfile << "$$LABEL\/1,part1" << std::endl << "$$DATE\/070920" << std::endl;
	outfile << "$$DIMENSION\/" << std::fixed << std::setprecision(6) << std::setw(8) << std::setfill('0') <<
		left[0] << "," << left[1] << "," << left[2] << "," << right[0] << "," << right[1] << "," << right[2] << endl;
	outfile << "$$LAYERS\/" << layer.size() << endl;
	outfile << "$$HEADEREND";

	// write geometry
	for (int i(0); i < layer.size(); ++i) {
		ContourEdge* l = layer[i];
		map<int, int> mp = mpvector[i];
		int* lpedgeNum = lpenumvector[i];
		double layerThick = layerHeights[i];
		//        if (i == 0) layerThick = layerHeights[i];
		//        else layerThick = layerHeights[i] - layerHeights[i - 1];
		outfile << "$$LAYER\/" << layerThick << std::endl;
		size_t dir, num;
		double area;
		double p[3];
		for (auto it = mp.begin(); it != mp.end(); it++) {
			int loopId = it->first; int eid = it->second;
			num = lpedgeNum[loopId];
			if (num < 3) continue;
			area = areaCal(l, eid);
			if (abs(area) < 1e-6) continue;
			dir = area > 0 ? 1 : 0;
			outfile << "$$POLYLINE\/1," << dir << "," << num + 1 << ",";
			//输出这个loop中的所有点
			ContourEdge fe = l[eid];
			p[0] = fe.stpt[0]; p[1] = fe.stpt[1]; p[2] = fe.stpt[2];
			outfile /*<< std::fixed*/ << p[0] << "," << p[1] << ",";
			int nextId = fe.next;
			while (nextId != eid) {
				ContourEdge ne = l[nextId];
				p[0] = ne.stpt[0]; p[1] = ne.stpt[1]; p[2] = ne.stpt[2];
				outfile /*<< std::fixed*/ << p[0] << "," << p[1] << ",";
				nextId = ne.next;
			}
			p[0] = fe.stpt[0]; p[1] = fe.stpt[1]; p[2] = fe.stpt[2];
			outfile /*<< std::fixed*/ << p[0] << "," << p[1] << std::endl;
		}

	}
	outfile.close();
	return true;
}
















