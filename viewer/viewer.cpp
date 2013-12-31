/****************************************************************************\

  viewer.c: an OpenGL model viewer
      adapted from AnyView by Cliff Woolley
  http://www.cs.virginia.edu/~jcc5t/classes/grad/animation/cart_pole/

\****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include "plylib.h"
#include "viewer.h"

#define FALSE 0
#define TRUE  1

/******************************************************************************
File I/O interface
******************************************************************************/

PlyProperty vert_props[] = { /* list of property information for a vertex */
  {"x",             PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,coord[X]),  0,0,0,0},
  {"y",             PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,coord[Y]),  0,0,0,0},
  {"z",             PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,coord[Z]),  0,0,0,0},
  {"nx",            PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,normal[X]), 0,0,0,0},
  {"ny",            PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,normal[Y]), 0,0,0,0},
  {"nz",            PLY_FLOAT, PLY_FLOAT, offsetof(Vertex,normal[Z]), 0,0,0,0},
};

PlyProperty face_props[] = { /* list of property information for a face */
  {"vertex_indices", PLY_INT, PLY_UINT, offsetof(Face,verts),
   1, PLY_UCHAR, PLY_UCHAR, offsetof(Face,nverts)},
};

#define verify_range(i,v)                                               \
        if ((v < 0) || (v >= (unsigned int)model->nverts)) {            \
            printf("Illegal vertex index %d at pos %d\n",v,i);          \
            exit(-1);                                                   \
        }

static void read_ply(FILE *fp, Model *model)
{
    int i,j;
    PlyFile *ply;
    int nprops;
    int num_elems;
    PlyProperty **plist;
    char *elem_name;
    int nelems;
    char **element_list;
    int has_vertices = FALSE, has_faces = FALSE;
    
    /*** Read in the original PLY object ***/
    ply  = ply_read(fp, &nelems, &element_list);

    for (i = 0; i < nelems; i++) {
        
        /* get the description of the first element */
        elem_name = element_list[i];
        plist = ply_get_element_description (ply, elem_name,
                                             &num_elems, &nprops);
        
        if (equal_strings ("vertex", elem_name)) {
            
            int has_x  = FALSE, has_y  = FALSE, has_z  = FALSE;
            int has_nx = FALSE, has_ny = FALSE, has_nz = FALSE;

            for (j=0; j<nprops; j++) {
                if (equal_strings("x", plist[j]->name)) {
                    ply_get_property (ply, elem_name, &vert_props[0]);  /* x */
                    has_x = TRUE;
                }
                else if (equal_strings("y", plist[j]->name)) {
                    ply_get_property (ply, elem_name, &vert_props[1]);  /* y */
                    has_y = TRUE;
                }
                else if (equal_strings("z", plist[j]->name)) {
                    ply_get_property (ply, elem_name, &vert_props[2]);  /* z */
                    has_z = TRUE;
                }
                else if (equal_strings("nx", plist[j]->name)) {
                    ply_get_property (ply, elem_name, &vert_props[3]);  /* nx*/
                    has_nx = TRUE;
                }
                else if (equal_strings("ny", plist[j]->name)) {
                    ply_get_property (ply, elem_name, &vert_props[4]);  /* ny*/
                    has_ny = TRUE;
                }
                else if (equal_strings("nz", plist[j]->name)) {
                    ply_get_property (ply, elem_name, &vert_props[5]);  /* nz*/
                    has_nz = TRUE;
                }
            }
            
            /* test for necessary properties */
            if (!(has_x && has_y && has_z)) {
                printf("Vertices don't have x, y, and z.\n");
                exit(-1);
            }
            has_vertices = TRUE;
            if (has_nx && has_ny && has_nz) {
                model->has_normals = TRUE;
            }

            /* create a vertex list to hold all the vertices */
            model->vertices = (Vertex *)malloc(num_elems*sizeof(Vertex));
            model->nverts = num_elems;
            
            /* grab all the vertex elements */
            for (j = 0; j < num_elems; j++) {
                ply_get_element (ply, &model->vertices[j]);
            }

        }
        
        else if (equal_strings ("face", elem_name)) {

            int k,l, breathing_room=64, alloced_faces, has_fverts = FALSE;
            Face face;

            for (j=0; j<nprops; j++) {
                if (equal_strings("vertex_indices", plist[j]->name)) {
                    /* vertex_indices */
                    ply_get_property(ply, elem_name, &face_props[0]);
                    has_fverts = TRUE;
                }
            }
            
            /* test for necessary properties */
            if (!has_fverts) {
                printf("Faces must have vertex indices.\n");
                exit(-1);
            }
            has_faces = TRUE;

            /* allocate enough space to handle all faces being triangles */
            model->fvlist =
                (unsigned int *)malloc(num_elems*3*sizeof(unsigned int));
            alloced_faces = num_elems;

            /* grab all the face elements */
            model->nfaces = 0;
            for (j=0,k=0; j < num_elems; j++) {
                ply_get_element (ply, &face);
                /* if we have any polygons with more than three
                 * verts, we're going to need more space. */
                model->nfaces += face.nverts - 2;
                if (model->nfaces > alloced_faces) {
                    breathing_room *= 2;
                    alloced_faces += breathing_room;
                    model->fvlist = (unsigned int *)realloc(model->fvlist,
                                         alloced_faces*3*sizeof(unsigned int));
                }
                /* this splits any polygons with more than
                 * three vertices into a set of triangles */
                verify_range(0,face.verts[0]);
                for (l=0; l<face.nverts-2; l++) {
                    model->fvlist[k++] = face.verts[0];
                    model->fvlist[k++] = face.verts[l+1];
                    model->fvlist[k++] = face.verts[l+2];
                    verify_range(l+1,face.verts[l+1]);
                    verify_range(l+2,face.verts[l+2]);
                }
                free(face.verts);
            }
        }

    }

    if (!(has_vertices && has_faces)) {
        printf("Missing required element(s).\n");
        exit(-1);
    }

    ply_close(ply);
}

void LoadModel(char *filename, Model *model)
{
    FILE *fp = fopen(filename, "rb");
    char buf[10];

    if (fp == NULL) {
        printf("Unable to open file.\n");
        exit(-1);
    }

    fgets(buf, sizeof(buf), fp);
    rewind(fp);

    if (strncmp(buf, "ply", 3) == 0) {
        read_ply(fp, model);
    }
}

void DrawModel(Model *model)
{
    if (model->has_normals) {
        glInterleavedArrays(GL_N3F_V3F, 0, model->vertices);
    }
    else {
        glDisableClientState(GL_NORMAL_ARRAY);
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, 6*sizeof(float),
                        (char *)model->vertices+3*sizeof(float));
    }
    glDrawElements(GL_TRIANGLES, model->nfaces*3,
                   GL_UNSIGNED_INT, model->fvlist);
}
