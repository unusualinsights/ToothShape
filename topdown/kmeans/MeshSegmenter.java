import PointVector;
import TriangleMesh;
import java.lang.*;
import java.util.*;
import java.util.Vector;

class TrianglePair
{
	public int face1, face2;

	public TrianglePair () {}
};

class MeshSegmenter
{
	TriangleMesh trimesh;
	int[][] neighbors;

	public MeshSegmenter (TriangleMesh trimesh)
	{
		this.trimesh = trimesh;

		// Compute edgeshare array
		int N = trimesh.nfaces;
		neighbors = new int[N][3];
		for (int i = 0; i < N; i++)
		{
			int curIndex = 0;
			for (int j = 0; j < N; j++)
			{
				if (i != j && ShareEdge (i, j))
				{
					neighbors[i][curIndex] = j;
					curIndex++;
				}
			}
		}
	}

	public Vector KmeansSegmentation (int k, int numits)
	{
		Hashtable segments = new Hashtable ();

		// facesIn[i] is the Vector of face numbers (Integers) lying in
		// cluster i
		Vector[] facesIn = new Vector[k];
		for (int i = 0; i < k; i++)
			facesIn[i] = new Vector ();

		// centers[i] is the center point of cluster i
		Point[] centers = new Point[k];
		for (int i = 0; i < k; i++)
			centers[i] = new Point ();

		// Pick k initial centers
		for (int i = 0; i < k; i++)
		{
			double f = Math.random ();
			int fnum = (int)(f * trimesh.nfaces);
			centers[i] = Centroid (trimesh.faces[fnum]);
		}

		// Cluster all faces into the corresponding k clusters
		for (int i = 0; i < trimesh.nfaces; i++)
		{
			double mindist = Double.MAX_VALUE;
			for (int j = 0; j < k; j++)
			{
				double dist = Point.Difference
				              (Centroid (trimesh.faces[i]),
				               centers[j]).Length ();
				if (dist < mindist)
				{
					mindist = dist;
					facesIn[j].addElement (new Integer (i));
				}
			}
		}

		// Repeatedly compute centers and recluster
		for (int i = 0; i < numits; i++)
		{
			// Compute the centroid of each of the k clusters
			for (int j = 0; j < k; j++)
			{
				Integer K0 = (Integer)facesIn[j].elementAt (0);
				int k0 = K0.intValue ();
				Point cen0 = Centroid (trimesh.faces[k0]);
				Point centroid = cen0;
				double sf = 1.0 / (double)facesIn[j].size ();
				for (int m = 1; m < facesIn[j].size (); m++)
				{
					int f = ((Integer)facesIn[j].elementAt
					                      (m)).intValue ();
					Point cen = Centroid (trimesh.faces[f]);
					Vec vec = Point.Difference (cen, cen0);
					vec.Scale (sf);
					centroid.Translate (vec);
				}
				centers[j] = centroid;
			}

			// Reinitialize clustering
			facesIn = new Vector[k];
			for (int j = 0; j < k; j++)
				facesIn[j] = new Vector ();

			// Compute corresponding clusters using new centers
			for (int j = 0; j < trimesh.nfaces; j++)
			{
				double mindist = Double.MAX_VALUE;
				for (int m = 0; m < k; m++)
				{
					double dist = Point.Difference
					       (Centroid (trimesh.faces[j]),
					        centers[m]).Length ();
					if (dist < mindist)
					{
						mindist = dist;
						facesIn[m].addElement
						           (new Integer (j));
					}
				}
			}
		}

		// Store clustering into Hashtable
		for (int cnum = 0; cnum < k; cnum++)
		{
			Integer Cnum = new Integer (cnum);
			for (int a = 0; a < facesIn[cnum].size (); a++)
			{
				Integer Fnum = (Integer)facesIn[cnum].elementAt
				                                      (a);
				segments.put (Fnum, Cnum);
			}
		}

		Vector segInfo = new Vector ();
		segInfo.addElement (segments);
		segInfo.addElement (new Integer (k));

		return segInfo;
	}

	// Euclidean distance between centroids of the two faces
	public double GetDistance (int face1, int face2)
	{
		Face tri1 = trimesh.faces[face1];
		Face tri2 = trimesh.faces[face2];
		Point c1 = Centroid (tri1);
		Point c2 = Centroid (tri2);
		return Point.Difference (c1, c2).Length ();
	}

	public double GetDistance (TrianglePair pair)
	{
		return GetDistance (pair.face1, pair.face2);
	}

	public Point Centroid (Face face)
	{
		Point p0 = trimesh.vertices[face.fvlist[0]].Coord ();
		Point p1 = trimesh.vertices[face.fvlist[1]].Coord ();
		Point p2 = trimesh.vertices[face.fvlist[2]].Coord ();
		Vec p1_p0 = Point.Difference (p1, p0);
		Vec p2_p0 = Point.Difference (p2, p0);
		p1_p0.Scale (1.0 / 3.0);
		p2_p0.Scale (1.0 / 3.0);
		return Point.Translate (p0, Vec.Add (p1_p0, p2_p0));
	}

	public TrianglePair[] GetPairs ()
	{
		int numPairs = 3 * trimesh.nfaces / 2;
		TrianglePair[] pairs = new TrianglePair[numPairs];

		int pairNum = 0;
		for (int i = 0; i < trimesh.nfaces; i++)
		{
			for (int j = i + 1; j < trimesh.nfaces; j++)
			{
				if (ShareEdge (i, j))
				{
					pairs[pairNum] = new TrianglePair ();
					pairs[pairNum].face1 = i;
					pairs[pairNum].face2 = j;
					pairNum++;
				}
			}
		}

		return pairs;
	}

	// Returns true iff faces[i] and faces[j] share an edge
	public boolean ShareEdge (int i, int j)
	{
		boolean first = false, second = false, third = false;
		first  = trimesh.faces[i].fvlist[0] ==
		         trimesh.faces[j].fvlist[0] ||
		         trimesh.faces[i].fvlist[0] ==
		         trimesh.faces[j].fvlist[1] ||
		         trimesh.faces[i].fvlist[0] ==
		         trimesh.faces[j].fvlist[2];
		second = trimesh.faces[i].fvlist[1] ==
		         trimesh.faces[j].fvlist[0] ||
		         trimesh.faces[i].fvlist[1] ==
		         trimesh.faces[j].fvlist[1] ||
		         trimesh.faces[i].fvlist[1] ==
		         trimesh.faces[j].fvlist[2];
		third  = trimesh.faces[i].fvlist[2] ==
		         trimesh.faces[j].fvlist[0] ||
		         trimesh.faces[i].fvlist[2] ==
		         trimesh.faces[j].fvlist[1] ||
		         trimesh.faces[i].fvlist[2] ==
		         trimesh.faces[j].fvlist[2];
		return (first && second || first && third || second && third);
	}

	public Vec UnnormalizedFaceNormal (Face triangle)
	{
		Vertex p0 = trimesh.vertices[triangle.fvlist[0]];
		Vertex p1 = trimesh.vertices[triangle.fvlist[1]];
		Vertex p2 = trimesh.vertices[triangle.fvlist[2]];
		Vec v1_0 = Point.Difference (p1.Coord (), p0.Coord ());
		Vec v2_0 = Point.Difference (p2.Coord (), p0.Coord ());
		return v1_0.Cross (v2_0);
	}

	// Returns [0], [1] as common vertices
	// [2] and [3] are the other vertices of the pair,
	// where [2] is in face 1 and [3] is in face 2 of the pair
	public Vertex[] Vertices (TrianglePair pair)
	{
		Face tri1 = trimesh.faces[pair.face1];
		Face tri2 = trimesh.faces[pair.face2];
		Vertex[] verts1 = new Vertex[3];
		Vertex[] verts2 = new Vertex[3];
		verts1[0] = trimesh.vertices[tri1.fvlist[0]];
		verts1[1] = trimesh.vertices[tri1.fvlist[1]];
		verts1[2] = trimesh.vertices[tri1.fvlist[2]];
		verts2[0] = trimesh.vertices[tri2.fvlist[0]];
		verts2[1] = trimesh.vertices[tri2.fvlist[1]];
		verts2[2] = trimesh.vertices[tri2.fvlist[2]];
		boolean[] equalVerts = new boolean[9];
		Vertex[] retverts = new Vertex[4];
		equalVerts[0] = (verts1[0] == verts2[0]);
		equalVerts[1] = (verts1[0] == verts2[1]);
		equalVerts[2] = (verts1[0] == verts2[2]);
		equalVerts[3] = (verts1[1] == verts2[0]);
		equalVerts[4] = (verts1[1] == verts2[1]);
		equalVerts[5] = (verts1[1] == verts2[2]);
		equalVerts[6] = (verts1[2] == verts2[0]);
		equalVerts[7] = (verts1[2] == verts2[1]);
		equalVerts[8] = (verts1[2] == verts2[2]);
		int eqto1in1 = 0, eqto1in2 = 0, eqto2in1 = 0, eqto2in2 = 0;
		int i;
		for (i = 0; i < 9; i++)
		{
			if (equalVerts[i])
			{
				eqto1in1 = i / 3;
				eqto1in2 = i % 3;
				break;
			}
		}
		for (i = 8; i >= 0; i--)
		{
			if (equalVerts[i])
			{
				eqto2in1 = i / 3;
				eqto2in2 = i % 3;
				break;
			}
		}
		retverts[0] = verts1[eqto1in1];
		retverts[1] = verts1[eqto2in1];
		int lone1, lone2;
		if (eqto1in1 == 0)
		{
			if (eqto2in1 == 1)
				lone1 = 2;
			else // eqto2in1 == 2
				lone1 = 1;
		}
		else if (eqto1in1 == 1)
		{
			if (eqto2in1 == 0)
				lone1 = 2;
			else // eqto2in1 == 2
				lone1 = 0;
		}
		else // eqto1in1 == 2
		{
			if (eqto2in1 == 0)
				lone1 = 1;
			else // eqto2in1 == 1
				lone1 = 0;
		}
		if (eqto1in2 == 0)
		{
			if (eqto2in2 == 1)
				lone2 = 2;
			else // eqto2in2 == 2
				lone2 = 1;
		}
		else if (eqto1in2 == 1)
		{
			if (eqto2in2 == 0)
				lone2 = 2;
			else // eqto2in2 == 2
				lone2 = 0;
		}
		else // eqto1in2 == 2
		{
			if (eqto2in2 == 0)
				lone2 = 1;
			else // eqto2in2 == 1
				lone2 = 0;
		}
		retverts[2] = verts1[lone1];
		retverts[3] = verts2[lone2];
		return retverts;
	}
};




