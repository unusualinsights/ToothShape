import PointVector;
import TriangleMesh;
import Matrix;
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

	// neighbors[i] is set of indices of faces neighboring trimesh.faces[i]
	int[][] neighbors;

	// starfaces[i] is set of indices of faces touching trimesh.vertices[i]
	int[][] starfaces;

	// k[i] is curvature at vertex i
	double[] k;

	// onormals[i] is normal vector at vertex i
	Vec[] onormals;

	public MeshSegmenter (TriangleMesh trimesh)
	{
		this.trimesh = trimesh;

		k = new double[trimesh.nverts];

		// Compute edgeshare array
		neighbors = new int[trimesh.nfaces][3];
		for (int i = 0; i < trimesh.nfaces; i++)
		{
			int curIndex = 0;
			for (int j = 0; j < trimesh.nfaces; j++)
			{
				if (i != j && ShareEdge (i, j))
				{
					neighbors[i][curIndex] = j;
					curIndex++;
				}
			}
		}

		Hashtable stars = new Hashtable ();
		for (int facenum = 0; facenum < trimesh.nfaces; facenum++)
		{
			int v0 = trimesh.faces[facenum].fvlist[0];
			int v1 = trimesh.faces[facenum].fvlist[1];
			int v2 = trimesh.faces[facenum].fvlist[2];

			Vector V0 = (Vector)stars.get (new Integer (v0));
			if (V0 == null) V0 = new Vector ();
			Vector V1 = (Vector)stars.get (new Integer (v1));
			if (V1 == null) V1 = new Vector ();
			Vector V2 = (Vector)stars.get (new Integer (v2));
			if (V2 == null) V2 = new Vector ();

			Integer Facenum = new Integer (facenum);
			V0.addElement (Facenum);
			V1.addElement (Facenum);
			V2.addElement (Facenum);

			stars.put (new Integer (v0), V0);
			stars.put (new Integer (v1), V1);
			stars.put (new Integer (v2), V2);
		}
		starfaces = new int[trimesh.nverts][];
		for (Enumeration e = stars.keys (); e.hasMoreElements ();)
		{
			Integer Vnum = (Integer)e.nextElement ();
			Vector Faces = (Vector)stars.get (Vnum);
			int vnum = Vnum.intValue ();
			starfaces[vnum] = new int[Faces.size ()];
			for (int i = 0; i < Faces.size (); i++)
			{
				Integer Fnum = (Integer)Faces.elementAt (i);
				starfaces[vnum][i] = Fnum.intValue ();
			}
		}

		onormals = new Vec[trimesh.nfaces];
		for (int i = 0; i < trimesh.nfaces; i++)
			onormals[i] = OrientedFaceNormal (i);
	}

	// See "Partitioning 3D Surface Meshes Using Watershed Segmentation"
	// by Alan P. Mangan and Ross T. Whitaker, IEEE Viz Vol 5(4), 1999
	public Vector WatershedSegmentation (double threshold)
	{
		Hashtable segments = new Hashtable ();

		// 1. Compute principal curvatures and directions
		for (int i = 0; i < trimesh.nverts; i++)
		{
			k[i] = CovarianceCurvature (i);
		}

		// 2. Label each local minimum
		boolean[] label = new boolean[trimesh.nverts];
		for (int i = 0; i < trimesh.nverts; i++)
			label[i] = IsLocalMinimum (i);

		// Convert vertex segmentation to face segmentation
		// ------------------------------------------------
		// If any vertex of triangle fnum is labeled, then that
		// triangle lies in segment zero, else in segment one.
		// A labeled vertex is one which lies on a "segment"
		// boundary.  But in this algorithm, we are treating every
		// triangle touching a "segment" boundary as *segment 0* and
		// all other triangles as *segment 1*.
		int numSegs = 2;
		for (int fnum = 0; fnum < trimesh.nfaces; fnum++)
		{
			boolean s0 = label[trimesh.faces[fnum].fvlist[0]];
			boolean s1 = label[trimesh.faces[fnum].fvlist[1]];
			boolean s2 = label[trimesh.faces[fnum].fvlist[2]];
			int snum = (s0 || s1 || s2) ? 0 : 1;
			segments.put (new Integer (fnum), new Integer (snum));
		}

		Vector segInfo = new Vector ();
		segInfo.addElement (segments);
		segInfo.addElement (new Integer (numSegs));

		return segInfo;
	}

	public double CovarianceCurvature (int vnum)
	{
		int N = starfaces[vnum].length; // # faces in star of vnum
		double oneOverN = 1.0 / N;
		double[][] C = new double[3][3]; // covariance matrix

		// Initialize all entries in C to 0
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				C[i][j] = 0.0;

		// Compute face normals and average normal
		Vec[] fnorms = new Vec[N];
		Vec avgNormal = new Vec ();
		for (int i = 0; i < N; i++)
		{
			fnorms[i] = onormals[starfaces[vnum][i]];
			avgNormal.Add (fnorms[i]);
		}
		avgNormal.Scale (oneOverN);
		double x_ = avgNormal.GetX ();
		double y_ = avgNormal.GetY ();
		double z_ = avgNormal.GetZ ();

		// Compute (co)variances
		for (int i = 0; i < N; i++)
		{
			double xt = fnorms[i].GetX ();
			double yt = fnorms[i].GetY ();
			double zt = fnorms[i].GetZ ();

			double xt_x_ = xt - x_;
			double yt_y_ = yt - y_;
			double zt_z_ = zt - z_;

			C[0][0] += xt_x_ * xt_x_;
			C[1][1] += yt_y_ * yt_y_;
			C[2][2] += zt_z_ * zt_z_;

			C[0][1] += xt_x_ * yt_y_;
			C[1][0] += yt_y_ * xt_x_;

			C[0][2] += xt_x_ * zt_z_;
			C[2][0] += zt_z_ * xt_x_;

			C[1][2] += yt_y_ * zt_z_;
			C[2][1] += zt_z_ * yt_y_;
		}
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				C[i][j] *= oneOverN;
				C[i][j] = Math.sqrt (C[i][j]);
			}

		// Curvature is norm of covariance matrix
		Matrix M = new Matrix (C);
		return M.normF ();
	}

	public boolean IsLocalMinimum (int vnum)
	{
		// Look at all neighbors of vnum, check whether k[vnum] is
		// strictly less than k of all neighbors of vnum
		Vector Neighbors = new Vector ();
		for (int i = 0; i < starfaces[vnum].length; i++)
		{
			int fnbr = starfaces[vnum][i];
			for (int j = 0; j < 3; j++)
			{
				int nbr = trimesh.faces[fnbr].fvlist[j];
				if (nbr == vnum) continue;
				Integer Nbr = new Integer (nbr);
				if (!Neighbors.contains (Nbr))
				{
					Neighbors.addElement (Nbr);
					if (!(k[vnum] < k[nbr]))
					// vnum is NOT a local minimum
						return false;
				}
			}
		}

		// vnum IS a local minimum
		return true;
	}

	public Vec NormalizedWeightedNormalAtVertex (int vnum)
	{
		Vec normal = new Vec ();

		// Get oriented normals of faces in starfaces[vnum]
		Vec[] facenormals = new Vec[starfaces[vnum].length];
		for (int i = 0; i < starfaces[vnum].length; i++)
			facenormals[i]=OrientedFaceNormal(starfaces[vnum][i]);

		// Set normal equal to normalized weighted sum of those normals
		for (int i = 0; i < starfaces[vnum].length; i++)
		{
			Vec normToAdd = facenormals[i];
			normToAdd.Scale (AreaOfFace (starfaces[vnum][i]));
			normal.Add (normToAdd);
		}
		normal.Normalize ();

		return normal;
	}

	public Vec FaceNormal (int fnum)
	{
		int v0 = trimesh.faces[fnum].fvlist[0];
		int v1 = trimesh.faces[fnum].fvlist[1];
		int v2 = trimesh.faces[fnum].fvlist[2];
		Point p0 = trimesh.vertices[v0].Coord ();
		Point p1 = trimesh.vertices[v1].Coord ();
		Point p2 = trimesh.vertices[v2].Coord ();
		Vec p1_p0 = Point.Difference (p1, p0);
		Vec p2_p0 = Point.Difference (p2, p0);
		Vec normal = p1_p0.Cross (p2_p0);
		normal.Normalize ();
		return normal;
	}

	public Vec OrientedFaceNormal (int fnum)
	{
		Vec normal = FaceNormal (fnum);
		int parity = NumIntersections (Centroid (fnum), normal);
		if (parity % 2 == 1)
			normal.Scale (-1.0);
		return normal;
	}

	public int NumIntersections (Point o, Vec d)
	{
		int numint = 0;
		for (int fnum = 0; fnum < trimesh.nfaces; fnum++)
			if (RayIntersects (fnum, o, d))
				numint++;
		return numint;
	}

	public boolean RayIntersects (int fnum, Point o, Vec d)
	{
		int v0 = trimesh.faces[fnum].fvlist[0];
		int v1 = trimesh.faces[fnum].fvlist[1];
		int v2 = trimesh.faces[fnum].fvlist[2];
		Point p0 = trimesh.vertices[v0].Coord ();
		Point p1 = trimesh.vertices[v1].Coord ();
		Point p2 = trimesh.vertices[v2].Coord ();
		Vec e1 = Point.Difference (p1, p0);
		Vec e2 = Point.Difference (p2, p0);
		Vec s1 = d.Cross (e2);
		double divisor = s1.Dot (e1);
		if (divisor == 0.0)
			return false;
		double invDivisor = 1.0 / divisor;

		// See if first barycentric coordinate lies in triangle
		Vec s = Point.Difference (o, p0);
		double b1 = s.Dot (s1) * invDivisor;
		if (b1 < 0.0 || b1 > 1.0)
			return false;

		// See if second barycentric coordinate lies in triangle
		Vec s2 = s.Cross (e1);
		double b2 = d.Dot (s2) * invDivisor;
		if (b2 < 0.0 || b1 + b2 > 1.0)
			return false;

		// Compute t to intersection point
		double t = e2.Dot (s2) * invDivisor;
		if (t <= 0)
			return false;

		return true;
	}

	public double AreaOfFace (int fnum)
	{
		int v0 = trimesh.faces[fnum].fvlist[0];
		int v1 = trimesh.faces[fnum].fvlist[1];
		int v2 = trimesh.faces[fnum].fvlist[2];
		Point p0 = trimesh.vertices[v0].Coord ();
		Point p1 = trimesh.vertices[v1].Coord ();
		Point p2 = trimesh.vertices[v2].Coord ();
		Vec p1_p0 = Point.Difference (p1, p0);
		Vec p2_p0 = Point.Difference (p2, p0);
		// Area of triangle (p0, p1, p2) is 1/2 ||p1_p0 x p2_p0||
		return 0.5 * Math.abs (p1_p0.Cross (p2_p0).Length ());
	}

	// Returns distances between centroids of *adjacent* triangles
	// These are the "edge weights" of the dual graph of the triangle mesh
	public double[] GetDistances ()
	{
		int numPairs = 3 * trimesh.nfaces / 2;
		double[] pairdist = new double[numPairs];

		int pairNum = 0;
		for (int i = 0; i < trimesh.nfaces; i++)
		{
			for (int j = i + 1; j < trimesh.nfaces; j++)
			{
				if (ShareEdge (i, j))
				{
					pairdist[pairNum] = GetDistance (i, j);
					pairNum++;
				}
			}
		}

		return pairdist;
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

	public Point Centroid (int fnum)
	{
		return Centroid (trimesh.faces[fnum]);
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




