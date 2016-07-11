#pragma once

#include "glm/glm.hpp"
#include "Types.h"
#include "Scene.h"
#include <ctime>
#include <iostream>

#include "string"
#include "atlimage.h"

#define G 6.674e-11
#define C 3e+8
#define PI 3.141592
#define Angle 45*PI/180
#define DirEps 0.01
#define II 256
#define JJ 225
#define RHO 4

class CTracer
{
public:
  SRay MakeRay(glm::uvec2 pixelPos, double xPoint, double yPoint);  // Create ray for specified pixel
  glm::vec3 TraceRay(SRay& ray, CImage *fon, CImage *texture); // Trace ray, compute its color
  void RenderImage(int xRes, int yRes);
  void SaveImageToFile(std::string fileName);
  CImage* LoadImageFromFile(std::string fileName);
  glm::vec3 Antialiacing(int j, int i, CImage *fon, CImage *texture);

public:
  SCamera m_camera;
  CScene* m_pScene;
  int Antial;
  int AlphaD;
  int Sph;
  CTracer(){
	Antial = 0;
	AlphaD = 0;
	Sph = 0;
  }
};

struct CPhoton{
	glm::vec3 m_dir;
	glm::vec3 m_pos;
	glm::vec3 m_acc;
	glm::vec3 m_dst;
	void print_photon() {
		printf("dir = %lf, %lf, %lf \n", m_dir.x, m_dir.y, m_dir.z);
		printf("pos = %lf, %lf, %lf \n", m_pos.x, m_pos.y, m_pos.z);
		printf("acc = %lf, %lf, %lf \n", m_acc.x, m_acc.y, m_acc.z);
		printf("dst = %lf, %lf, %lf \n\n", m_dst.x, m_dst.y, m_dst.z);
	}
};