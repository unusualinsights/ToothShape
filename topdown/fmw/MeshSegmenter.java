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

	// Info for curvature tensor for each vertex
	double[] k1;
	double[] k2;
	Vec   [] T1;
	Vec   [] T2;

	// Normal vector at each vertex
	Vec[] Normals;

	public MeshSegmenter (TriangleMesh trimesh)
	{
		this.trimesh = trimesh;

		k1 = new double[trimesh.nverts];
		k2 = new double[trimesh.nverts];
		T1 = new Vec   [trimesh.nverts];
		T2 = new Vec   [trimesh.nverts];

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

		Normals = new Vec[trimesh.nverts];
	}

	// See "Perception-based 3D Triangle Mesh Segmentation Using Fast
	// Marching Watersheds" by D. L. Page, A. F. Koschan, and M. A. Abidi
	// http://imaging.utk.edu/~page/MyResearch/pubs/page-cvpr2003-158.pdf
	// Author suggestion: use alpha = 0.3 for most applications
	//                    and nhdRadius = 1
	public Vector FastMarchingWatersheds (double alpha, int nhdRadius)
	{
		Hashtable segments = new Hashtable ();

		// 1. Compute principal curvatures and directions
		for (int i = 0; i < trimesh.nverts; i++)
		{
			Hashtable tensor = GetCurvatureTensor (i);
			k1[i] = ((Double)tensor.get ("k1")).doubleValue ();
			k2[i] = ((Double)tensor.get ("k2")).doubleValue ();
			T1[i] = (Vec)tensor.get ("T1");
			T2[i] = (Vec)tensor.get ("T2");
			Normals[i] = (Vec)tensor.get ("Normal");
//System.out.println ("Vertex #" + i + ":");
//System.out.println ("k1[" + i + "] = " + k1[i]);
//System.out.println ("k2[" + i + "] = " + k2[i]);
//System.out.println ("T1[" + i + "] = " + T1[i]);
//System.out.println ("T2[" + i + "] = " + T2[i]);
		}

		// 2. Threshold each minimum k2 in mesh
		double thresh = 0.0;
		Vector F = new Vector (); // indices of feature vertices
		int numNegVerts = 0;
		for (int i = 0; i < trimesh.nverts; i++)
		{
			if (k2[i] <= 0.0)
			{
//System.out.println ("k2[" + i + "] = " + k2[i]);
				thresh += k2[i];
				numNegVerts++;
			}
		}
		thresh *= alpha / (double)numNegVerts;
//System.out.println ("threshold = " + thresh);

		for (int i = 0; i < trimesh.nverts; i++)
			if (k2[i] <= thresh)
				F.addElement (new Integer (i));

		// 3. Apply mathematical morphology to clean up segmentation
		//    Result is marker set F_
		Vector F_ = Open (Close (F, nhdRadius), nhdRadius);
		//Vector F_ = F;
		boolean[] label = new boolean[trimesh.nverts];
		for (int i = 0; i < trimesh.nverts; i++)
			label[i] = F_.contains (new Integer (i));

		// 4. Grow marker set to cover entire mesh
/*
		MaxHeap heap = new MaxHeap ();
		for (int i = 0; i < trimesh.nverts; i++)
			if (label[i])
				heap = Grow (i, heap, label);
		while (!heap.isEmpty ())
		{
			FMWData data = heap.Pop ();
			int i = data.vertex;
			if (!label[i])
			{
				label[i] = data.label;
				heap = Grow (i, heap, label);
			}
		}
*/

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

	public MaxHeap Grow (int vnum, MaxHeap heap, boolean[] label)
	{
		Vector nbrs = Nhd (vnum);
		for (int i = 0; i < nbrs.size (); i++)
		{
			int j = ((Integer)nbrs.elementAt (i)).intValue ();
			if (!label[j])
			{
				double key = DirectionalHeight (vnum, j);
				FMWData data = new FMWData (j, label[vnum]);
				heap.Push (key, data);
			}
		}
		return heap;
	}

	public double DirectionalHeight (int i, int j)
	{
		Matrix I = Matrix.identity (3, 3);
		Matrix Nj = Normals[j].toColMatrix ();
		Matrix NjNjt = Nj.times (Nj.transpose ());
		Matrix Tj = I.minus (NjNjt);
		Point vi = trimesh.vertices[i].Coord ();
		Point vj = trimesh.vertices[j].Coord ();
		Vec vi_vj = Point.Difference (vi, vj);
		Matrix TJI = Tj.times (vi_vj.toColMatrix ());
		TJI.timesEquals (1.0 / vi_vj.Length ());
		Vec Tji = new Vec (TJI.get (0, 0),
		                   TJI.get (1, 0),
		                   TJI.get (2, 0));
		double cos2thetai = Tji.Dot (T1[j]);
		cos2thetai *= cos2thetai;
		double sin2thetai = 1.0 - cos2thetai;
		double kji = k1[j] * cos2thetai + k2[j] * sin2thetai;
		return kji;
	}

	public Vector Nhd (int vnum)
	{
		Vector nbrs = new Vector ();
		for (int i = 0; i < starfaces[vnum].length; i++)
		{
			int fnum = starfaces[vnum][i];
			Face face = trimesh.faces[fnum];
			Integer V0 = new Integer (face.fvlist[0]);
			Integer V1 = new Integer (face.fvlist[1]);
			Integer V2 = new Integer (face.fvlist[2]);
			if (!nbrs.contains (V0)) nbrs.addElement (V0);
			if (!nbrs.contains (V1)) nbrs.addElement (V1);
			if (!nbrs.contains (V2)) nbrs.addElement (V2);
		}
		return nbrs;
	}

	public Vector Nhd (Vector F)
	{
		Vector G = F;
		for (int i = 0; i < F.size (); i++)
		{
			int v = ((Integer)F.elementAt (i)).intValue ();
			Vector nbrs = Nhd (v);
			for (int j = 0; j < nbrs.size (); j++)
			{
				Integer nbr = (Integer)nbrs.elementAt (j);
				if (!G.contains (nbr)) G.addElement (nbr);
			}
		}
		return G;
	}

	public Vector Nhd (Vector F, int n)
	{
		if (n < 1) return null;
		if (n == 1) return Nhd (F);
		return Nhd (Nhd (F, n - 1));
	}

	public Vector Dilate (Vector F, int n)
	{
		return Nhd (F, n);
	}

	public Vector Erode (Vector F, int n)
	{
		if (n < 1) return null;
		if (n > 1) return Erode (Erode (F, n - 1), 1);

		Vector G = new Vector ();
		for (int i = 0; i < F.size (); i++)
		{
			Integer V = (Integer)F.elementAt (i);
			Vector vnbrs = Nhd (V.intValue ());
			boolean allInF = true;
			for (int j = 0; j < vnbrs.size (); j++)
			{
				Integer Nbr = (Integer)vnbrs.elementAt (j);
				if (!F.contains (Nbr))
				{
					allInF = false;
					break;
				}
			}
			if (allInF) G.addElement (V);
		}
		return G;
	}

	public Vector Open (Vector F, int n)
	{
		return Dilate (Erode (F, n), n);
	}

	public Vector Close (Vector F, int n)
	{
		return Erode (Dilate (F, n), n);
	}

	public Hashtable GetCurvatureTensor (int vnum)
	{
		// 1a. Estimate normal vector at vertex
		Vec normal = NormalizedWeightedNormalAtVertex (vnum);
//System.out.println ("Nv" + vnum + ": " + normal.toString ());

		// 1b. Estimate curvature tensor matrix

		// Compute neighboring vertices' indices
		Vector Neighbors = new Vector ();
		for (int i = 0; i < starfaces[vnum].length; i++)
		{
			int facenum = starfaces[vnum][i];
			int v0 = trimesh.faces[facenum].fvlist[0];
			int v1 = trimesh.faces[facenum].fvlist[1];
			int v2 = trimesh.faces[facenum].fvlist[2];
			Integer V0 = new Integer (v0);
			Integer V1 = new Integer (v1);
			Integer V2 = new Integer (v2);
			if (v0 != vnum && !Neighbors.contains (V0))
				Neighbors.addElement (V0);
			if (v1 != vnum && !Neighbors.contains (V1))
				Neighbors.addElement (V1);
			if (v2 != vnum && !Neighbors.contains (V2))
				Neighbors.addElement (V2);
		}

		// Compute weights of neighboring vertices
		Matrix M = new Matrix (3, 3);
		double[] w = new double[Neighbors.size ()];
		double sum = 0.0;
		for (int i = 0; i < Neighbors.size (); i++)
		{
			int j = ((Integer)Neighbors.elementAt (i)).intValue ();
			int first = 0, second = 0;
			for (int k = 0; k < starfaces[vnum].length; k++)
			{
				for (int m = 0; m < 3; m++)
				{
					int f = starfaces[vnum][k];
					if (trimesh.faces[f].fvlist[m] == j)
						first = starfaces[vnum][k];
				}
			}
			for (int k = starfaces[vnum].length - 1; k >= 0; k--)
			{
				for (int m = 0; m < 3; m++)
				{
					int f = starfaces[vnum][k];
					if (trimesh.faces[f].fvlist[m] == j)
						second = starfaces[vnum][k];
				}
			}
			w[i] = AreaOfFace (first) + AreaOfFace (second);
			sum += w[i];
		}
		double sumrecip = 1.0 / sum;
		for (int i = 0; i < Neighbors.size (); i++)
			w[i] *= sumrecip;
//Print out weights
//for (int i = 0; i < Neighbors.size (); i++)
//{
//int j = ((Integer)Neighbors.elementAt (i)).intValue ();
//System.out.println ("w" + vnum + j + " = " + w[i]);
//}

		// Compute the matrix
		Matrix I = Matrix.identity (3, 3);
		for (int i = 0; i < Neighbors.size (); i++)
		{
			int j = ((Integer)Neighbors.elementAt (i)).intValue ();
			Matrix Nv = normal.toColMatrix ();
			Matrix Nvt = normal.toRowMatrix ();
			Matrix dif = I.minus (Nv.times (Nvt));
			Point vvnum = trimesh.vertices[vnum].Coord ();
			Point vj = trimesh.vertices[j].Coord ();
			Vec vec = Point.Difference (vj, vvnum);
			Matrix vv = vec.toColMatrix ();
			Matrix T = dif.times (vv);
//System.out.println ("T[" + vnum + "][" + j + "] = (" + T.get (0, 0) + ", " +
//T.get (1, 0) + ", " + T.get (2, 0) + ")");
//System.out.println ("TNorm" + vnum + j + " = " + T.normF ());
			double norm = T.normF ();
			if (norm != 0.0)
				T.timesEquals (1.0 / norm);
			Matrix prod = T.times (T.transpose ());
			Vec vj_vi = Point.Difference (vj, vvnum);
			double k = 2.0 * normal.Dot (vj_vi);
			k /= vj_vi.LengthSquared ();
//System.out.println ("k[" + vnum + "][" + j + "] = " + k);
			prod.timesEquals (w[i] * k);
//if (vnum == 0 && j == 1)
//{
//System.out.println ("w01k01T01T01t = ");
//for (int aa = 0; aa < 3; aa++)
//{
//for (int bb = 0; bb < 3; bb++)
//System.out.print (prod.get (aa, bb) + " ");
//System.out.println ();
//}
//}
			M.plusEquals (prod);
		}
//System.out.println ("Matrix M_v" + vnum + ":");
//for (int i = 0; i < 3; i++)
//for (int j = 0; j < 3; j++)
//System.out.println (M.get (i, j));

		// 1c. Obtain final eigenvalues and eigenvectors

		// Compute matrix for Householder transformation
		Vec E1 = new Vec (1.0, 0.0, 0.0);
		Vec WW = new Vec ();
		Vec plus = Vec.Add (E1, normal);
		Vec minus = Vec.Subtract (E1, normal);
		if (minus.LengthSquared () > plus.LengthSquared ())
			WW = minus;
		else
			WW = plus;
		WW.Normalize ();
		Matrix W = WW.toColMatrix ();

		// Householder matrix
		Matrix Q = I.minus (W.times (W.transpose ()).times (2.0));
//System.out.println ("Matrix Q_v" + vnum + ":");
//for (int i = 0; i < 3; i++)
//for (int j = 0; j < 3; j++)
//System.out.println (Q.get (i, j));

		// Orthonormal basis vectors for tangent space at vertex
		Vec T1_ = new Vec (Q.get (0, 1), Q.get (1, 1), Q.get (2, 1));
		Vec T2_ = new Vec (Q.get (0, 2), Q.get (1, 2), Q.get (2, 2));

		Matrix QtMQ = Q.transpose ().times (M).times (Q);
		// Check that first row and column of QtMQ are zero
		// Check that QtMQ[2][1] == QtMQ[1][2]
//System.out.println ("Matrix Q_v" + vnum + "^tM_v" + vnum + "Q_v" + vnum+":");
//for (int i = 0; i < 3; i++)
//for (int j = 0; j < 3; j++)
//System.out.println (QtMQ.get (i, j));
//Matrix MQ = M.times (Q);
//System.out.println ("Matrix M_v" + vnum + "Q_v" + vnum+":");
//for (int i = 0; i < 3; i++)
//for (int j = 0; j < 3; j++)
//System.out.println (MQ.get (i, j));
		double m11_ = QtMQ.get (1, 1);
		double m12_ = QtMQ.get (1, 2); // equals QtMQ.get (2, 1)
		double m22_ = QtMQ.get (2, 2);

//System.out.println ();
//System.out.println ("Vertex " + vnum + ":");
		// Diagonalize QtMQ with Givens rotation in tangent plane
		double x = normal.GetX ();
		double y = normal.GetY ();
		double z = normal.GetZ ();
		double m11m12_m22_ = m11_ * m12_ / m22_;
		double a = y * z * m11_ + (z * z - 1) * m12_ +
                           (y * y - 1) * m11m12_m22_ + y * z * m11_;
		double b = m12_ + m11m12_m22_;
		double m11 = m11_, m22 = m22_;
		Vec T1 = T1_;
		Vec T2 = T2_;
		// Theorem: a == 0 iff b == 0
		// If a == 0 and b == 0, then the matrix is already
		// diagonalized, i.e. m11_ = m11 and m22_ == m22.  Otherwise
		// need to diagonalize matrix with a Givens rotation:
		if (a != 0)
		{
			double costheta = 1.0 + b / a;
			double sintheta = Math.sqrt (1 - costheta * costheta);
			m11 = (1.0 + (1.0 - costheta) * (y * y - 1)) * m11_ +
			      (-x * sintheta + (1.0 - costheta) * y * z) * m12_;
			m22 = (x * sintheta + (1.0 - costheta) * y * z) * m12_ +
			      (1.0 + (1.0 - costheta) * (z * z - 1)) * m22_;
			T1 = Vec.Subtract (Vec.Scale (T1_, costheta),
			                       Vec.Scale (T2_, sintheta));
			T2 = Vec.Add      (Vec.Scale (T1_, sintheta),
			                       Vec.Scale (T2_, costheta));
		}
		Double k1 = new Double (3.0 * m11 - m22);
		Double k2 = new Double (3.0 * m22 - m11);
//System.out.println ("m11 = " + m11);
//System.out.println ("m22 = " + m22);

		// Return (k1, k2, T1, T2)
		Hashtable H = new Hashtable ();
		H.put ("k1", k1);
		H.put ("k2", k2);
		H.put ("T1", T1);
		H.put ("T2", T2);
		H.put ("Normal", normal);
		return H;
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




