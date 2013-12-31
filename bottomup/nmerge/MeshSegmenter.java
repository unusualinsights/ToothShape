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

class SegmentPair
{
	public int seg1, seg2;

	public SegmentPair () {}
};

class MeshSegmenter
{
	TriangleMesh trimesh;

	public MeshSegmenter (TriangleMesh trimesh)
	{
		this.trimesh = trimesh;
	}

	// Merges pairs up to and including thresholdScore
	public Vector PairMergeSegmentation (double thresholdScore)
	{
		Hashtable segments = new Hashtable ();

		// Initially each triangle is its own segment
		for (int i = 0; i < trimesh.nfaces; i++)
		{
			Integer I = new Integer (i);
			segments.put (I, I);
		}

		// Compute scores of adjacent pairs of triangles
		int numPairs = 3 * trimesh.nfaces / 2; // # pairs = # edges
		TrianglePair[] pairs = GetPairs ();
		Hashtable mergeScores = new Hashtable ();
		for (int i = 0; i < numPairs; i++)
		{
			SegmentPair pair = new SegmentPair ();
			pair.seg1 = pairs[i].face1;
			pair.seg2 = pairs[i].face2;
			Double score = new Double (NormalScore (pairs[i]));
			mergeScores.put (pair, score);
		}

		// Merge segments up to threshold score
		int numSegments = trimesh.nfaces;
		while (true)
		{
			// Get the pair to merge with minimum score
			SegmentPair pair = ExtractMin (mergeScores);
			Double Score = (Double)mergeScores.get (pair);

			// Check if we are done segmenting
			if (Score.doubleValue () < thresholdScore ||
			    mergeScores.size () < 1) break;

			// Remove pair (one, two) from mergeScores and replace
			// occurrences of two with one, avoiding duplicates
			mergeScores = RemoveAndUpdate
			              (pair, mergeScores, segments);

			// Merge segments pair.seg1 and pair.seg2
			segments = UpdateSegments (pair, segments);
			numSegments--;
		}

		Vector segInfo = new Vector ();
		segInfo.addElement (segments);
		segInfo.addElement (new Integer (numSegments));
		segInfo.addElement (new Integer (trimesh.nfaces));

		return segInfo;
	}

	public SegmentPair ExtractMin (Hashtable table)
	{
		SegmentPair minpair = new SegmentPair ();
		minpair.seg1 = -1;
		minpair.seg2 = -1;
		double min = MAX_DOUBLE;
		for (Enumeration e = table.keys (); e.hasMoreElements ();)
		{
			SegmentPair p = (SegmentPair)e.nextElement ();
			Double score = (Double)table.get (p);
			double s = score.doubleValue ();
			if (s < min)
			{
				min = s;
				minpair = p;
			}
		}

		return minpair;
	}

	// If pair = (one, two), then all triangles in segment two become
	// triangles in segment one and the pair (one, two) is removed; must
	// also remove any duplicates in the Hashtable, though that is kinda
	// not allowed by Java to have two entries with same key in Hashtable
	public Hashtable RemoveAndUpdate (SegmentPair pair,
	                                  Hashtable mergeScores,
	                                  Hashtable segments)
	{
		Hashtable table = mergeScores;
		table.remove (pair);
		int one = pair.seg1;
		int two = pair.seg2;
		for (Enumeration e = table.keys (); e.hasMoreElements ();)
		{
			// Get next pair p from table
			SegmentPair p = (SegmentPair)e.nextElement ();

			// Change any twos to ones in pair p
			if (p.seg1 == two)
				p.seg1 = one;
			if (p.seg2 == two)
				p.seg2 = one;

			// If p used to contain two and now contains one in
			// its place, or if p originally contained one, must
			// recompute its score since one and two together form
			// the new segment we will call one
			if (p.seg1 == one || p.seg2 == one)
			{
				p = GetPair (p, table);
				table.remove (p);
				// Compute updated score of merging pair p, but
				// only if p is a valid pair--in particular, p
				// should not be of the form (x, x), since it
				// does not make sense to merge a segment with
				// itself
				if (p.seg1 != p.seg2)
				{
					Double s = new Double (NmergeScore
					                       (p, segments));
					table.put (p, s);
				}
			}
		}

		return table;
	}

	public static boolean TableContains (SegmentPair pair, Hashtable table)
	{
		SegmentPair reverse = ReversePair (pair);
		return table.contains (pair) || table.contains (reverse);
	}

	public static SegmentPair GetPair (SegmentPair pair, Hashtable table)
	{
		SegmentPair reverse = ReversePair (pair);
		if (table.contains (pair))
			return pair;
		if (!table.contains (reverse))
			return null;
		return reverse;
	}

	public static SegmentPair ReversePair (SegmentPair pair)
	{
		SegmentPair newpair = new SegmentPair ();
		newpair.seg1 = pair.seg2;
		newpair.seg2 = pair.seg1;
		return newpair;
	}

	public Hashtable UpdateSegments (SegmentPair pair, Hashtable segments)
	{
		Hashtable h = segments;
		Integer one = new Integer (pair.seg1);
		Integer two = new Integer (pair.seg2);

		// Replace all entries of the form (x, two) with (x, one)
		for (Enumeration e = h.keys (); e.hasMoreElements ();)
		{
			Integer trinum = (Integer)e.nextElement ();
			Integer segnum = (Integer)h.get (trinum);
			if (segnum.equals (two))
			{
				h.remove (trinum);
				h.put (trinum, one);
			}
		}

		return h;
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
		first  = trimesh.faces[i].fvlist[0]
		         == trimesh.faces[j].fvlist[0] ||
		         trimesh.faces[i].fvlist[0]
		         == trimesh.faces[j].fvlist[1] ||
		         trimesh.faces[i].fvlist[0]
		         == trimesh.faces[j].fvlist[2];
		second = trimesh.faces[i].fvlist[1]
		         == trimesh.faces[j].fvlist[0] ||
		         trimesh.faces[i].fvlist[1]
		         == trimesh.faces[j].fvlist[1] ||
		         trimesh.faces[i].fvlist[1]
		         == trimesh.faces[j].fvlist[2];
		third  = trimesh.faces[i].fvlist[2]
		         == trimesh.faces[j].fvlist[0] ||
		         trimesh.faces[i].fvlist[2]
		         == trimesh.faces[j].fvlist[1] ||
		         trimesh.faces[i].fvlist[2]
		         == trimesh.faces[j].fvlist[2];
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

	public double NormalScore (TrianglePair pair)
	{
		Vec[] normals = new Vec[2];
		Vertex[] p = Vertices (pair);
		Vec p1_p0 = Point.Difference (p[1].Coord (), p[0].Coord ());
		Vec p2_p0 = Point.Difference (p[2].Coord (), p[0].Coord ());
		Vec p3_p0 = Point.Difference (p[3].Coord (), p[0].Coord ());
		normals[0] = p1_p0.Cross (p2_p0);
		normals[1] = p3_p0.Cross (p1_p0);
		normals[0].Normalize ();
		normals[1].Normalize ();
		return normals[0].Dot (normals[1]);
	}

	public double NmergeScore (SegmentPair pair, Hashtable segments)
	{
		Vector One = new Vector (); // indices of triangles in seg. 1
		Vector Two = new Vector (); // indices of triangles in seg. 2
		Integer one = new Integer (pair.seg1);
		Integer two = new Integer (pair.seg2);

		// Get indices of triangles in segments one and two
		for (Enumeration e = segments.keys (); e.hasMoreElements ();)
		{
			Integer trinum = (Integer)e.nextElement ();
			Integer segnum = (Integer)segments.get (trinum);
			if (segnum.equals (one)) One.addElement (trinum);
			if (segnum.equals (two)) Two.addElement (trinum);
		}

		return NmergeScore (One, Two);
	}

	public double NmergeScore (Vector One, Vector Two)
	{
		Vector I = Union (One, Two);
		Vector[] normals = new Vector[I.size ()];

		// Find first face and normal
		int i0 = (int)(Math.random () * I.size ());
		int ind0 = ((Integer)I.get (i0)).intValue ();
		normals[0] = OrientedFaceNormal (ind0);
		Face face0 = trimesh.faces[ind0];
		Point c0 = Centroid (face0);
		I.removeElement (new Integer (ind0));

		Integer[] Indices = new Integer[I.size ()];
		I.copyInto (Indices);
		double[] cdist = new double[Indices.length];
		// Compute array of centroid distances to reference centroid
		for (int i = 0; i < Indices.length; i++)
		{
			int indi = Indices[i].intValue ();
			Point ci = Centroid (trimesh.faces[indi]);
			cdist[i] = Point.Difference (ci, c0).Length ();
		}

		// Sort centroid distances and corresponding indices
		NormalQsort (cdist, Indices);

		// Compute remaining faces/normals
		for (int i = 1; i < Indices.length + 1; i++)
		{
			int indi = Indices[i].intValue ();
			normals[i] = OrientedFaceNormal (indi);
		}

		// Measure amount of self-overlap in the array of normals
		return Overlap (normals, I);
	}

	// Amount of self-overlap in ordered Vector of normals
	// I is Vector of indices, where
	// normals[i] is normal vector of triangle number
	//     (Integer)(I.elementAt (i)).intValue ()
	public double Overlap (Vector[] normals, Vector I)
	{
		// First 3 normals form first triangle on sphere
		for (int i = 3; i < normals.length; i++)
		{
			if (ith face normal intersects area covered)
			{
				overlap += amount of overlap in steradians
			}
		}
		return overlap;
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

	// Returns Vector of A union B
	public static Vector Union (Vector A, Vector B)
	{
		Vector C = A;
		for (int i = 0; i < B.size (); i++)
			A.addElement (B.elementAt (i));
		return C;
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

	// Sort cdist and corresponding Indices in descending order of cdist,
	// i.e. postcondition implies
	// cdist[0] >= cdist[1] >= ... >= cdist[cdist.length - 1]
	void NormalQsort (double[] cdist, Integer[] Indices)
	{
		NormalQsort (cdist, Indices, 0, cdist.length - 1);
	}

	void NormalQsort (double[] qq, Integer[] pp, int low, int high)
	{
		int lo = low;
		int hi = high;
		if (lo >= hi) return;
		double mid = qq[(lo + hi) / 2];
		while (lo <= hi)
		{
			while (lo < high && mid < qq[lo])
				lo++;
			while (hi > low && mid > qq[hi])
				hi--;
			if (lo <= hi)
			{
				double w = qq[lo];
				qq[lo] = qq[hi];
				qq[hi] = w;
				Integer ind = pp[lo];
				pp[lo] = pp[hi];
				pp[hi] = ind;
				lo++;
				hi--;
			}
		}

		if (low < hi)
			NormalQsort (qq, pp, low, hi);
		if (lo < high)
			NormalQsort (qq, pp. lo, high);
	}
};




