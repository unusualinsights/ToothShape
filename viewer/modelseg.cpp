
#include <assert.h>
#include "modelseg.h"

//
// Returns uniform random number in [0, 1]
//
float ranf ()
{
        double  UniformRandomNum;
                                                                                
        UniformRandomNum = (double)rand();
        UniformRandomNum /= (double)RAND_MAX;
        return (float)UniformRandomNum; // from 0 to 1
}

//
// Function definitions for class ModelSegmenter
//

ModelSegmenter::ModelSegmenter ()
{
}

ModelSegmenter::ModelSegmenter (const Model &m)
{
	int i, N = m.nverts;
	model.has_normals = m.has_normals;
	model.nverts = N;
	model.nfaces = m.nfaces;
	//model.vertices = (Vertex *)malloc (N * sizeof (Vertex));
	for (i = 0; i < N; i++)
		model.vertices[i] = m.vertices[i];
	//model.fvlist = (unsigned int *)malloc (N * 3 * sizeof (unsigned int));
	for (i = 0; i < 3 * m.nfaces; i++)
		model.fvlist[i] = m.fvlist[i];
}

void ModelSegmenter::NaiveSegmentation
(
	int fvlist[MAX_NUMBER_OF_SEGMENTS][MAX_NUMBER_OF_TRIANGLES * 3 + 1],
	int &numSegments
)
const
{
	// Randomly pick a vertex
	// Collect every triangle in its star into one segment
	// Repeat until every triangle has been added, or:
	//   numSegments == MAX_NUMBER_OF_SEGMENTS - 1
	//   in which case we add all remaining triangles to one last segment

	int i;
	int numVerts = model.nverts, numTris = model.nfaces;
	const int TSIZE = numTris * 3;
	const int FVSIZE = MAX_NUMBER_OF_TRIANGLES * 3 + 1;
	Vertex **V = new Vertex*[numVerts];
	int *T = new int[TSIZE];
	for (i = 0; i < numVerts; i++)
	{
		V[i] = new Vertex;
		*V[i] = model.vertices[i];
	}
	for (i = 0; i < numTris; i++)
	{
		assert (0 <= 3 * i);
		assert (3 * i + 2 < TSIZE);
		assert (0 <= model.fvlist[3 * i]);
		assert (model.fvlist[3 * i] < numVerts);
		assert (0 <= model.fvlist[3 * i + 1]);
		assert (model.fvlist[3 * i + 1] < numVerts);
		assert (0 <= model.fvlist[3 * i + 2]);
		assert (model.fvlist[3 * i + 2] < numVerts);
		T[  3 * i  ] = model.fvlist[  3 * i  ];
		T[3 * i + 1] = model.fvlist[3 * i + 1];
		T[3 * i + 2] = model.fvlist[3 * i + 2];
	}

	int numTrisAdded = 0;

	for (i = 0; i < MAX_NUMBER_OF_SEGMENTS; i++)
	{
		if (i < MAX_NUMBER_OF_SEGMENTS - 1)
		{
			int curIndex = 0;
			int k;
			for (k = 0; k < numVerts; k++) if (V[k] != NULL) break;
			assert (k < numVerts);
			// Fill in fvlist[i][] with triangles in star of V[k]
			for (int j = 0; j < numTris; j++)
			{
				assert (0 <= 3 * j);
				assert (3 * j + 2 < TSIZE);
				// T[3j..3j+2] touches V[k] -> add to fvlist[i]
				if
				(
					k == T[  3 * j  ] ||
					k == T[3 * j + 1] ||
					k == T[3 * j + 2]
				)
				{
					// Add T[3j..3j+2] to fvlist[i][]
					assert (0 <= i);
					assert (i < MAX_NUMBER_OF_SEGMENTS);
					assert (0 <= curIndex);
					assert (curIndex + 2 < FVSIZE);
					fvlist[i][curIndex++] = T[  3 * j  ];
					fvlist[i][curIndex++] = T[3 * j + 1];
					fvlist[i][curIndex++] = T[3 * j + 2];
					assert (0 <= T[3 * j]);
					assert (T[3 * j] < numVerts);
					V[T[  3 * j  ]] = NULL;
					assert (0 <= T[3 * j + 1]);
					assert (T[3 * j + 1] < numVerts);
					V[T[3 * j + 1]] = NULL;
					assert (0 <= T[3 * j + 2]);
					assert (T[3 * j + 2] < numVerts);
					V[T[3 * j + 2]] = NULL;
					T[  3 * j  ] = -1;
					T[3 * j + 1] = -1;
					T[3 * j + 2] = -1;
					numTrisAdded++;
				}
			}
			// 3 * this number is the # of ints in fvlist[i][]
			fvlist[i][MAX_NUMBER_OF_TRIANGLES * 3] = numTrisAdded;
			// If all triangles have been added, break
			if (numTrisAdded == numTris) break;
		}
		else
		{
			// Add remaining triangles in T to last segment
			int numTrisAdded = 0;
			int curIndex = 0;
			for (int j = 0; j < numTris; j++)
			{
				assert (0 <= 3 * j);
				assert (3 * j + 2 < TSIZE);
				if (T[3 * j] != -1)
				{
					assert (0 <= i);
					assert (i < MAX_NUMBER_OF_SEGMENTS);
					assert (0 <= curIndex);
					assert (T[3 * j + 1] != -1);
					assert (T[3 * j + 2] != -1);
					assert (curIndex + 2 < FVSIZE);
					fvlist[i][curIndex++] = T[  3 * j  ];
					fvlist[i][curIndex++] = T[3 * j + 1];
					fvlist[i][curIndex++] = T[3 * j + 2];
					assert (0 <= T[3 * j]);
if (T[3*j] >= numVerts) printf ("j = %d, T[3*j] = %d\n", j, T[3*j]);
					assert (T[3 * j] < numVerts);
					V[T[  3 * j  ]] = NULL;
					assert (0 <= T[3 * j + 1]);
					assert (T[3 * j + 1] < numVerts);
					V[T[3 * j + 1]] = NULL;
					assert (0 <= T[3 * j + 2]);
					assert (T[3 * j + 2] < numVerts);
					V[T[3 * j + 2]] = NULL;
					T[  3 * j  ] = -1;
					T[3 * j + 1] = -1;
					T[3 * j + 2] = -1;
					numTrisAdded++;
				}
			}
			// 3 * this number is the # of ints in fvlist[i][]
			// The number itself is the # of triangles in segment
			assert (0 <= i);
			assert (i < MAX_NUMBER_OF_SEGMENTS);
			fvlist[i][MAX_NUMBER_OF_TRIANGLES * 3] = numTrisAdded;
		}
	}
	numSegments = i;

	delete [] V;
	delete [] T;
}

void ModelSegmenter::ReadSegmentation
(
	char filename[256],
	int fvlist[MAX_NUMBER_OF_SEGMENTS][MAX_NUMBER_OF_TRIANGLES * 3 + 1],
	int &numSegments
)
const
{
	FILE *file;
	file = fopen (filename, "r");
	if (file == NULL)
	{
		printf ("Error: can't open file.\n");
		return;
	}

	fscanf (file, "%d", &numSegments);
	int i, j;
	for (i = 0; i < numSegments; i++)
	{
		int numtris;
		fscanf (file, "%d", &numtris);
		fvlist[i][MAX_NUMBER_OF_TRIANGLES * 3] = numtris;
		for (j = 0; j < numtris; j++)
		{
			int trinum;
			int fsf = fscanf (file, "%d", &trinum);
			assert (fsf > 0);
			fvlist[i][  3 * j  ] = model.fvlist[  3 * trinum  ];
			fvlist[i][3 * j + 1] = model.fvlist[3 * trinum + 1];
			fvlist[i][3 * j + 2] = model.fvlist[3 * trinum + 2];
		}
	}

	fclose (file);
}


