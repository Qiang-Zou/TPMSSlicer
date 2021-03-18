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
#include <experimental/unordered_map>
#include <time.h>

#include <QImage>


TPMSModeler::TPMSModeler() :
    tolerence(1e-6),
    pathOfFiles("../results/")
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
    case 1: ptrModel = std::static_pointer_cast<TPMS>(std::make_shared<TPMSSchwarz>(up, low, left, right));
        break;
    case 2: ptrModel = std::static_pointer_cast<TPMS>(std::make_shared<TPMSMultiscale>(up, low, left, right));
        break;
    case 3: ptrModel = std::static_pointer_cast<TPMS>(std::make_shared<TPMSMultiscale2>(up, low, left, right));
        break;
    case 4: ptrModel = std::static_pointer_cast<TPMS>(std::make_shared<TPMSMultiscaleDistFieldBound>(up, low, left, right, pathOfFiles + "bone.txt"));
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

    //    double layerThick = 10;
    //    size_t resX(800), resY(800);
    Eigen::Vector3d left, right;
    ptrModel->getBBox(left, right);

    long time = clock();

    // a functor for slicing at each iteration
    // the image stores in-pixel as 1 (only these pixels are stored) and out-pixel as 0
    auto generateBinImage = [&](double layerLevel, int layerNum, string path){
        // padding of the image to handle the openness of TPMS at boundaries
        Eigen::Vector3d leftLocal(left), rightLocal(right);
        double sizeX = right[0] - left[0], sizeY = right[1] - left[1];
        leftLocal += Eigen::Vector3d(-sizeX / 20, -sizeY / 20, 0);
        rightLocal += Eigen::Vector3d(sizeX / 20, sizeY / 20, 0);
        leftLocal[2] = rightLocal[2] = layerLevel;

        // determine the cell size
        double cellSize(std::min((rightLocal[0] - leftLocal[0]) / std::min(resX, resY), (rightLocal[1] - leftLocal[1]) / std::min(resX, resY)));
        size_t xRes(std::ceil((rightLocal[0] - leftLocal[0]) / cellSize));
        size_t yRes(std::ceil((rightLocal[1] - leftLocal[1]) / cellSize));
        if (xRes < yRes) yRes = std::max(yRes, std::max(resX, resY));
        else xRes = std::max(xRes, std::max(resX, resY));
        std::cout << "TPMSModeler > digital image resoltuion at "<< layerNum << "-th layer :" << xRes << "*" << yRes << endl;
        std::unordered_map<std::string, double> digitImage; // string: "rowIdx+space+colIdx"

        // compute each pixel's binary value
        for (int rowIdx(0); rowIdx <= xRes; ++rowIdx) {
            for (int colIdx(0); colIdx <= yRes; ++colIdx) {
                double x = leftLocal[0] + rowIdx * cellSize;
                double y = leftLocal[1] + colIdx * cellSize;
                double z = layerLevel;
                if (ptrModel->isInSolid(x, y, z)) {
                    auto key = std::to_string(rowIdx) + " " + std::to_string(colIdx);
                    digitImage.insert(std::make_pair(key, 1));
                }
            }
        }
        if(digitImage.empty()) return;

        // save the bindary image to file
        int width = xRes, height = yRes;
        QImage img(width, height, QImage::Format_RGB32);
        for (size_t i(0); i < width; ++i) {
            for (size_t j(0); j < height; ++j) {
                auto iter = digitImage.find(std::to_string(i) + " " + std::to_string(j));
                if (iter != digitImage.end()) img.setPixel(i, j, qRgb(0, 0, 0));
                else img.setPixel(i, j, qRgb(255, 255, 255));
            }
        }
        string path_ = pathOfFiles + "layer_" + std::to_string(layerNum) + ".jpg";
        if (!img.save(QString::fromStdString(path_), "JPEG", 50)) std::cerr << "TPMSModeler > error writing picture at layer: " << layerNum << endl;
        else std::cout << "TPMSModeler > successfully writing picture at layer: " << layerNum << endl;

        // do contouring with optimization
        double leftCorner[2] = {0, 0};
        MarchingSquareBin msb(digitImage, xRes, yRes, cellSize, layerLevel - left[2], leftCorner); // binImage will be moved, not copy, into MarchingSquareBin
        msb.doContouringWithOptim();
        //msb.doContouring();
        //path = pathFiles + "layer_" + std::to_string(layerNum) + ".off";
        //msb.writeContours(path);
        auto result = msb.getContours(false);
        //.......................
        // do things with "result", e.g., save to "path"

        //........................
        std::cout << "LatticeModeler > successfully generating contours at layer: " << layerNum << endl;
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



















