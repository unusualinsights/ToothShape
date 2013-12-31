/****************************************************************************\

  viewer.h: adapted from AnyView by Cliff Woolley
  http://www.cs.virginia.edu/~jcc5t/classes/grad/animation/cart_pole/

\****************************************************************************/

#define X 0
#define Y 1
#define Z 2

typedef float Vec3[3];
typedef float Point[3];

typedef struct Vertex {
    Vec3  normal;
    Point coord;
} Vertex;

typedef struct Face {
    unsigned int *verts;     /* vertex index list */
    unsigned char nverts;    /* number of vertex indices in list */
} Face;

typedef struct Model {
    int            nverts,nfaces;
    Vertex        *vertices; /* vertex data */
    unsigned int  *fvlist;   /* triangle vertex index array */
    int            has_normals;
} Model;

void LoadModel(char *filename, Model *model);
void DrawModel(Model *model);
