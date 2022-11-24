#pragma once

#include <queue>
#include <vector>
#include <Eigen/Sparse>

#include "../PolyMesh/include/PolyMesh/PolyMesh.h"

namespace freeb
{
	// return the edge is or not the cutting edge
	std::vector<char> CutEdge(acamcad::polymesh::PolyMesh* mesh);
}