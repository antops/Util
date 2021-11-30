/*
*	created by JinMing 2019/11/06
*   射线簇 文件IO接口
*   version
*/
//All in Global Geometry
#ifndef RAYS_H
#define RAYS_H
#include <vector>
#include "vector3.h"


class RayVec
{
public:
	RayVec() {};
	~RayVec() {};
	void InitialRays() {
		int RayNum = 50;
		Rays.resize(RayNum);
		for (int i = 0; i < RayNum; i++) {
			Rays[i].resize(2);
			for (int j = 0; j < 1; j++) {
				Rays[i][j].x = 0.0;			Rays[i][j].y = 0.0;							Rays[i][j].z = 0.0;
				Rays[i][j + 1].x = 0.0;		Rays[i][j + 1].y = -0.5 + 0.02*(i + 0.5);	Rays[i][j + 1].z = 1.0;
			}
		}
	}

	void AddPointToRay(int index, Vector3 point) {
		Rays[index].push_back(point);
	}

	void insertRay2RayVec(std::vector<Vector3> ray) {
		if (ray.size() > 0) {//拒绝加空射线
			Rays.push_back(ray);
		}
	}
	errno_t writeRays2BiFile(std::string filename) {
		FILE* file;
		errno_t err;
		err = fopen_s(&file, filename.c_str(), "wb");
		if (err != 0) return err;
		int RayNum = Rays.size();
		int PointNum;
		fwrite(&RayNum, sizeof(int), 1, file);
		if (RayNum > 0) {//如果没有光线就不用写内容了
			for (int i = 0; i < RayNum; i++) {
				PointNum = Rays[i].size();
				fwrite(&PointNum, sizeof(int), 1, file);
				for (int j = 0; j < PointNum; j++) {
					fwrite(&Rays[i][j].x, sizeof(double), 1, file);
					fwrite(&Rays[i][j].y, sizeof(double), 1, file);
					fwrite(&Rays[i][j].z, sizeof(double), 1, file);
				}
			}
		}
		fclose(file);
		return err;
	}
	errno_t readRaysfromBiFile(std::string filename) {
		FILE* file;
		errno_t err;
		std::vector<std::vector<Vector3>> tempRays;
		err = fopen_s(&file, filename.c_str(), "rb");
		if (err != 0) return err;
		int RayNum;
		fread(&RayNum, sizeof(int), 1, file);
		tempRays.resize(RayNum);
		if (RayNum > 0) {//不止有一条光线，申请内存都内容
			for (int i = 0; i < RayNum; i++) {
				int PointNum;
				Vector3 Point;
				fread(&PointNum, sizeof(int), 1, file);
				for (int j = 0; j < PointNum; j++) {
					fread(&Point.x, sizeof(double), 1, file);
					fread(&Point.y, sizeof(double), 1, file);
					fread(&Point.z, sizeof(double), 1, file);
					tempRays[i].push_back(Point);
				}
			}
			Rays = tempRays;
		}
		fclose(file);
		return err;
	}

	std::vector<std::vector<Vector3>> Rays;

};


#endif // !DRAW_H