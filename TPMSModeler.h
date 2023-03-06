#pragma once

#include <list>
#include <string>
#include <memory>
#include <deque>
#include <set>
#include <fstream>
#include <iostream>
#include <map>

#include "../GLKLib/GLKObList.h"

#include "QMeshLib/QMesh/QMeshNode.h"
#include "QMeshLib/QMesh/QMeshEdge.h"
#include "QMeshLib/QMesh/QMeshFace.h"
#include "QMeshLib/QMesh/QMeshPatch.h"

#include "Eigen/Core"
#include "ParallelMarchingSquare.h"
using namespace Eigen;

class TPMS;

class TPMSModeler
{
public:
    TPMSModeler();
    ~TPMSModeler();

    // TPMS geometry and file
    void readTPMSParameters();
    void readTPMSParameters(QMeshPatch *&ptrMesh);
    QMeshPatch* doMarchingCube();


    // 3d printing
    bool sliceModelImplicit();
	void outputOFF(string &path, int pixelNum, int* lpedgeNum, int loopNum, map<int, int> &mp, ContourEdge* edges);
	bool writeCLI(std::string &path, std::vector<ContourEdge*> &layer, vector<int> &edgeNum, vector<map<int, int>> &mpvector, vector<int*> &lpenumvector);

private:
    std::shared_ptr<TPMS> ptrModel;
    double tolerence;
    std::string pathOfFiles;
    std::vector<std::shared_ptr<QMeshPatch>> layers;
	std::vector<ContourEdge*> layer;
	std::vector<int> edgeNum;
	std::vector<map<int, int>>mpvector;
	std::vector<int*>lpenumvector;
};


///////////////////////////////////////////////////////////////////////////////
/// classes of tpms
///
///
///
///

const double PI = 3.1415926;

class TPMS {
public:
    TPMS() {}
    TPMS(double upperBound_, double lowerBound_) : upperBound(upperBound_), lowerBound(lowerBound_) {}
    TPMS(int idx_,double upperBound_, double lowerBound_, Vector3d left_, Vector3d right_) : idx(idx_),upperBound(upperBound_), lowerBound(lowerBound_), left(left_), right(right_) {}
    ~TPMS(){}

    virtual std::string getSurfaceType() = 0;
    virtual double getValue(double x, double y, double z) = 0;
    virtual bool isInSolid(double x, double y, double z) {
        if (!isInBoundary(x, y, z)) return false;
        double lv = getValue(x, y, z);
        return lv <= upperBound && lv >= lowerBound;
    }

    virtual bool isInBoundary(double x, double y, double z) { return !(x < left[0] || x > right[0] || y < left[1] || y > right[1] || z < left[2] || z > right[2]);}
    void setBounds(double upperBound_, double lowerBound_) {upperBound = upperBound_; lowerBound = lowerBound_;}
    void getBounds(double& upperBound_, double& lowerBound_) {upperBound_ = upperBound; lowerBound_ = lowerBound;}
    void setBBox(Vector3d left_, Vector3d right_) {left = left_; right = right_;}
    void getBBox(Vector3d& left_, Vector3d& right_) {left_ = left; right_ = right;}
	void getIdx(int& idx_) { idx_ = idx; }

	int idx;
    double upperBound; // upperBound and lowerBound define the solid
    double lowerBound;
    Vector3d left; // bounding box
    Vector3d right;
};

class TPMSSchwarz : public TPMS{
public:
    TPMSSchwarz(double upperBound_, double lowerBound_) : TPMS(upperBound_, lowerBound_) {}
    TPMSSchwarz(int idx_,double upperBound_, double lowerBound_, Vector3d left_, Vector3d right_) : TPMS(idx_,upperBound_, lowerBound_, left_, right_) {}

    std::string getSurfaceType() {return std::string("Schwarz-Primitive");}

    double getValue(double x, double y, double z) {
        return cos(x) + cos(y) + cos(z);
    }
};

class TPMSMultiscale : public TPMS{
public:
    TPMSMultiscale(double upperBound_, double lowerBound_) : TPMS(upperBound_, lowerBound_) {}
    TPMSMultiscale(int idx_,double upperBound_, double lowerBound_, Vector3d left_, Vector3d right_) : TPMS(idx_,upperBound_, lowerBound_, left_, right_) {}

    std::string getSurfaceType() {return std::string("multiscale");}

    double getValue(double x, double y, double z) {
        double l1(10), l2(3);
        double w1(2 * PI / l1), w2(2 * PI / l2);
        auto phiG = cos(w1 * x) * sin(w1 * y) + cos(w1 * y) * sin(w1 * z) + cos(w1 * z) * sin(w1 * x);
        auto phiGSon = cos(w2 * x) * sin(w2 * y) + cos(w2 * y) * sin(w2 * z) + cos(w2 * z) * sin(w2 * x);
        return min(phiG, phiGSon);
    }
};

class TPMSMultiscale2 : public TPMS {
public:
    TPMSMultiscale2(double upperBound_, double lowerBound_) : TPMS(upperBound_, lowerBound_) {}
    TPMSMultiscale2(int idx_,double upperBound_, double lowerBound_, Vector3d left_, Vector3d right_) : TPMS(idx_,upperBound_, lowerBound_, left_, right_) {}

    std::string getSurfaceType() {return std::string("multiscale2");}

    // inner solid -> those less than 0
    double getValue(double x, double y, double z) {
        double l1(24), l2(6), l3(1.5);
        double w1(2 * PI / l1), w2(2 * PI / l2), w3(2 * PI / l3);
        auto Phi_G = cos(w1*x)*sin(w1*y)+cos(w1*y)*sin(w1*z)+cos(w1*z)*sin(w1*x)+0.7;
        auto Phi_G2 = cos(w1*x)*sin(w1*y)+cos(w1*y)*sin(w1*z)+cos(w1*z)*sin(w1*x)-0.7;
        auto Phi_G_son= cos(w2*x)*sin(w2*y)+cos(w2*y)*sin(w2*z)+cos(w2*z)*sin(w2*x)+0.7;
        auto Phi_G_son_2= cos(w2*x)*sin(w2*y)+cos(w2*y)*sin(w2*z)+cos(w2*z)*sin(w2*x)-0.7;
        auto Phi_G_son_3= cos(w3*x)*sin(w3*y)+cos(w3*y)*sin(w3*z)+cos(w3*z)*sin(w3*x)+0.5;
        auto Phi_G_son_3_2= cos(w3*x)*sin(w3*y)+cos(w3*y)*sin(w3*z)+cos(w3*z)*sin(w3*x)-0.5;
        // first and second level;
        auto Phi1=min(Phi_G,-Phi_G2);
        auto Phi2=min(Phi_G_son,-Phi_G_son_2);
        auto Phi_f=min(Phi1,Phi2);
        // second and third level;
        auto Phi3=min(Phi_G_son,-Phi_G_son_2);
        auto Phi4=min(Phi_G_son_3,-Phi_G_son_3_2);
        auto Phi_son=min(Phi3,Phi4);
        return min(Phi_f,Phi_son);
    }
};

class TPMSMultiscaleDistFieldBound : public TPMS{
public:
    TPMSMultiscaleDistFieldBound(double upperBound_, double lowerBound_, string path) : TPMS(upperBound_, lowerBound_) {setDistField(path);}
    TPMSMultiscaleDistFieldBound(int idx_,double upperBound_, double lowerBound_, Vector3d left_, Vector3d right_, string path) : TPMS(idx_,upperBound_, lowerBound_, left_, right_) {setDistField(path);}

    void setDistField(std::string path) {
        // open file
        std::ifstream ifs(path.c_str(), ifstream::in);
        std::string line;
        float height, xMin, yMin, xMax, yMax;
        int Nx, Ny;
        while (std::getline(ifs, line)) {
            stringstream str(line);
            str >> height >> xMin >> xMax >> yMin >> yMax >> Nx >> Ny;
            zHeight.push_back(height);
            bottomLeft.push_back(Vector2f(xMin, yMin)); upRight.push_back(Vector2f(xMax, yMax));
            nXY.push_back(make_pair(Nx + 1, Ny + 1)); // Nx and Ny are # of cells, grid has one more size then them
            resX.push_back((xMax - xMin) / Nx); resY.push_back((yMax - yMin) / Ny);
            // get distance values
            Nx++; Ny++;
            std::vector<float> entries(Nx * Ny, 0);
            for (auto& e : entries) str >> e;
            distField.emplace_back(entries);
        }
        ifs.close();
        resZ = (zHeight.back() - zHeight.front()) / (zHeight.size() - 1);
    }

    std::string getSurfaceType() {return std::string("TPMSMultiscaleDistFieldBound");}

    bool isInBoundary(double x, double y, double z) override {
        if (x < left[0] || x > right[0] || y < left[1] || y > right[1] || z < left[2] || z > right[2]) return false;
        // evaluate the distance sign
        if (z > zHeight.back() || z < zHeight.front()) return false;
        float zCoord = (z - zHeight.front()) / resZ;
        int zIdxLow = std::floor(zCoord);
        zCoord -= zIdxLow;

        if (x < bottomLeft[zIdxLow][0] || x > upRight[zIdxLow][0]) return false;
        float xCoord = (x - bottomLeft[zIdxLow][0]) / resX[zIdxLow];
        int xIdxLow = std::floor(xCoord);
        xCoord -= xIdxLow;

        if (y < bottomLeft[zIdxLow][1] || y > upRight[zIdxLow][1]) return false;
        float yCoord = (y - bottomLeft[zIdxLow][1]) / resY[zIdxLow];
        int yIdxLow = std::floor(yCoord);
        yCoord -= yIdxLow;

        int xIdxUP = (xIdxLow == nXY[zIdxLow].first)  ? xIdxLow : xIdxLow + 1;
        int yIdxUP = (yIdxLow == nXY[zIdxLow].second) ? yIdxLow : yIdxLow + 1;
        int zIdxUP = (zIdxLow == zHeight.size() - 1) ? zIdxLow : zIdxLow + 1;

        // Trilinear interpolation
        float c000 = distField[zIdxLow][xIdxLow + yIdxLow * nXY[zIdxLow].first];
        float c100 = distField[zIdxLow][xIdxUP + yIdxLow * nXY[zIdxLow].first];
        float c010 = distField[zIdxLow][xIdxLow + yIdxUP * nXY[zIdxLow].first];
        float c110 = distField[zIdxLow][xIdxUP  + yIdxUP * nXY[zIdxLow].first];
        float c001 = distField[zIdxUP][xIdxLow + yIdxLow * nXY[zIdxUP].first];
        float c101 = distField[zIdxUP][xIdxUP + yIdxLow  * nXY[zIdxUP].first];
        float c011 = distField[zIdxUP][xIdxLow + yIdxUP * nXY[zIdxUP].first];
        float c111 = distField[zIdxUP][xIdxUP + yIdxUP  * nXY[zIdxUP].first];

        float c00 = c000 * (1 - xCoord) + c100 * xCoord;
        float c01 = c001 * (1 - xCoord) + c101 * xCoord;
        float c10 = c010 * (1 - xCoord) + c110 * xCoord;
        float c11 = c011 * (1 - xCoord) + c111 * xCoord;

        float c0 = c00 * (1 - yCoord) + c10 * yCoord;
        float c1 = c01 * (1 - yCoord) + c11 * yCoord;

        float c = c0 * (1 - zCoord) + c1 * zCoord;
        return c <= 0;
    }

    bool isInSolid(double x, double y, double z) override {
        if (!isInBoundary(x, y, z)) return false;
        double lv = getValue(x, y, z);
        return lv <= upperBound && lv >= lowerBound;
    }


    // inner solid -> those less than 0
    double getValue(double x, double y, double z) {
        double l1(24), l2(6), l3(1.5);
        double w1(2 * PI / l1), w2(2 * PI / l2), w3(2 * PI / l3);
        auto Phi_G = cos(w1*x)*sin(w1*y)+cos(w1*y)*sin(w1*z)+cos(w1*z)*sin(w1*x)+0.7;
        auto Phi_G2 = cos(w1*x)*sin(w1*y)+cos(w1*y)*sin(w1*z)+cos(w1*z)*sin(w1*x)-0.7;
        auto Phi_G_son= cos(w2*x)*sin(w2*y)+cos(w2*y)*sin(w2*z)+cos(w2*z)*sin(w2*x)+0.7;
        auto Phi_G_son_2= cos(w2*x)*sin(w2*y)+cos(w2*y)*sin(w2*z)+cos(w2*z)*sin(w2*x)-0.7;
        auto Phi_G_son_3= cos(w3*x)*sin(w3*y)+cos(w3*y)*sin(w3*z)+cos(w3*z)*sin(w3*x)+0.5;
        auto Phi_G_son_3_2= cos(w3*x)*sin(w3*y)+cos(w3*y)*sin(w3*z)+cos(w3*z)*sin(w3*x)-0.5;
        // first and second level;
        auto Phi1=min(Phi_G,-Phi_G2);
        auto Phi2=min(Phi_G_son,-Phi_G_son_2);
        auto Phi_f=min(Phi1,Phi2);
        // second and third level;
        auto Phi3=min(Phi_G_son,-Phi_G_son_2);
        auto Phi4=min(Phi_G_son_3,-Phi_G_son_3_2);
        auto Phi_son=min(Phi3,Phi4);
        return min(Phi_f,Phi_son);
    }

private:
    std::vector<std::vector<float>> distField;
    std::vector<Vector2f> bottomLeft, upRight;
    std::vector<pair<int, int>> nXY;
    std::vector<float> resX, resY;
    std::vector<float> zHeight;
    float resZ;
};



