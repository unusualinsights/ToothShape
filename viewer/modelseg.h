//
// Specifications for classes for segmenting Models
// Models are the objects defined in viewer.h
//

#include <stdlib.h>
#include <stdio.h>
#include "viewer.h"

const int MAX_NUMBER_OF_SEGMENTS = 600;
const int MAX_NUMBER_OF_TRIANGLES = 6000;
const int MAX_NUMBER_OF_VERTICES = 5000;

//
// Model data type but with static array sizes
//
typedef struct StaticModel
{
	int nverts, nfaces;
	Vertex vertices[MAX_NUMBER_OF_VERTICES];
	unsigned int fvlist[MAX_NUMBER_OF_TRIANGLES * 3];
	int has_normals;
} StaticModel;

//
// Class for computing segmentations of a given Model
//
class ModelSegmenter
{
public:

	ModelSegmenter ();

	ModelSegmenter (const Model &m);

	void NaiveSegmentation
	(
		int fvlist[MAX_NUMBER_OF_SEGMENTS]
		          [MAX_NUMBER_OF_TRIANGLES * 3 + 1],
		int &numSegments
	)
	const;

	void ReadSegmentation
	(
		char filename[256],
		int fvlist[MAX_NUMBER_OF_SEGMENTS]
		          [MAX_NUMBER_OF_TRIANGLES * 3 + 1],
		int &numSegments
	)
	const;

private:

	StaticModel model;
};


