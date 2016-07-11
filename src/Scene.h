#pragma once

#include "Types.h"

struct CBlackHole{
	double m_mas;
	double m_rad;
	double m_radDisk;
	double m_phi;
	glm::vec3 m_pos;
};

class CScene
{
public:
	CBlackHole m_hole;
  // Set of meshes
};