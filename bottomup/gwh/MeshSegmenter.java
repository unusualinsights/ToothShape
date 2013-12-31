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

	public MeshSegmenter (TriangleMesh trimesh)
	{
		this.trimesh = trimesh;
	}

	// Merges pairs up to and including thresholdScore
	public Vector PairMergeSegmentation (double thresholdScore)
	{
		Hashtable segments = new Hashtable ();
		int i;

		// Initially each triangle is its own segment
		for (i = 0; i < trimesh.nfaces; i++)
		{
			Integer I = new Integer (i);
			segments.put (I, I);
		}

		// Compute scores of adjacent pairs of triangles
		int numPairs = 3 * trimesh.nfaces / 2; // # pairs = # edges
		TrianglePair[] pairs = GetPairs ();
		double[] scores = new double[numPairs];
		for (i = 0; i < numPairs; i++)
		{
			//scores[i] = i;
			scores[i] = NormalScore (pairs[i]);
		}

		// Sort triangle pair scores and pairs correspondingly
		PairScoreQuicksort (scores, pairs);

		// Check if list was actually sorted in descending order
		for (i = 0; i < numPairs - 1; i++)
		{
			if (scores[i] < scores[i + 1])
			{
				System.out.println ("Error: Scores not sorted");
				System.exit (0);
			}
		}

		int numSegments = trimesh.nfaces;
		// Merge segments up to threshold score
		for (i = 0; i < numPairs; i++)
		{
			if (scores[i] < thresholdScore) break;

			TrianglePair pair = pairs[i];

			// Check if pair.face1, face2 in same segment already
			Integer one = (Integer)segments.get
			              (new Integer (pair.face1));
			Integer two = (Integer)segments.get
			              (new Integer (pair.face2));
			if (one.equals (two)) continue;

			// Merge segments containing pair.face1, face2
			numSegments--;
			Integer newSegNum = one;
			for (Enumeration e = segments.keys ();
			     e.hasMoreElements ();)
			{
				Integer triNum = (Integer)e.nextElement ();
				Integer segNum = (Integer)segments.get (triNum);
				if (segNum.equals (two))
				{
					segments.remove (triNum);
					segments.put (triNum, newSegNum);
				}
			}
		}

		Vector segInfo = new Vector ();
		segInfo.addElement (segments);
		segInfo.addElement (new Integer (numSegments));
		segInfo.addElement (new Integer (trimesh.nfaces));

		return segInfo;
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
		first  = trimesh.faces[i].fvlist[0] == trimesh.faces[j].fvlist[0] ||
		         trimesh.faces[i].fvlist[0] == trimesh.faces[j].fvlist[1] ||
		         trimesh.faces[i].fvlist[0] == trimesh.faces[j].fvlist[2];
		second = trimesh.faces[i].fvlist[1] == trimesh.faces[j].fvlist[0] ||
		         trimesh.faces[i].fvlist[1] == trimesh.faces[j].fvlist[1] ||
		         trimesh.faces[i].fvlist[1] == trimesh.faces[j].fvlist[2];
		third  = trimesh.faces[i].fvlist[2] == trimesh.faces[j].fvlist[0] ||
		         trimesh.faces[i].fvlist[2] == trimesh.faces[j].fvlist[1] ||
		         trimesh.faces[i].fvlist[2] == trimesh.faces[j].fvlist[2];
		return (first && second || first && third || second && third);
	}

	public double AppaCurvatureScore (TrianglePair pair)
	{
		Face tri1 = trimesh.faces[pair.face1];
		Face tri2 = trimesh.faces[pair.face2];
/*
		Vertex p10 = trimesh.vertices[tri1.fvlist[0]];
		Vertex p11 = trimesh.vertices[tri1.fvlist[1]];
		Vertex p12 = trimesh.vertices[tri1.fvlist[2]];
		Vertex p20 = trimesh.vertices[tri2.fvlist[0]];
		Vertex p21 = trimesh.vertices[tri2.fvlist[1]];
		Vertex p22 = trimesh.vertices[tri2.fvlist[2]];
		Vec f1_1_0 = Point.Difference (p11.Coord (), p10.Coord ());
		Vec f1_2_0 = Point.Difference (p12.Coord (), p10.Coord ());
		Vec f2_1_0 = Point.Difference (p21.Coord (), p20.Coord ());
		Vec f2_2_0 = Point.Difference (p22.Coord (), p20.Coord ());
*/
		//f1_1_0.Normalize ();
		//f1_2_0.Normalize ();
		//f2_1_0.Normalize ();
		//f2_2_0.Normalize ();
		//Vec norm1 = f1_1_0.Cross (f1_2_0);
		//Vec norm2 = f2_1_0.Cross (f2_2_0);
		Vec norm1 = UnnormalizedFaceNormal (tri1);
		Vec norm2 = UnnormalizedFaceNormal (tri2);
		norm1.Normalize ();
		norm2.Normalize ();
		double dotproduct = norm1.Dot (norm2);
		return Math.abs (dotproduct);
		//return dotproduct;
		//return 1.0 - Math.abs (dotproduct);
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

	public double DumbUnusableCurvatureScore (TrianglePair pair)
	{
		// Set p1, p2 to be the common vertices
		// p3 and p4 the other vertices
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
		Vertex p1, p2, p3, p4;
		boolean[] equalVerts = new boolean[6];
		equalVerts[0] = (verts1[0] == verts2[0]);
		equalVerts[1] = (verts1[0] == verts2[1]);
		equalVerts[2] = (verts1[0] == verts2[2]);
		equalVerts[3] = (verts1[1] == verts2[1]);
		equalVerts[4] = (verts1[1] == verts2[2]);
		equalVerts[5] = (verts1[2] == verts2[2]);
		int eqto1in1 = 0, eqto1in2 = 0, eqto2in1 = 0, eqto2in2 = 0;
		int i;
		for (i = 0; i < 6; i++)
		{
			if (equalVerts[i])
			{
				eqto1in1 = (i < 3) ? 0 : ((i < 5) ? 1 : 2);
				eqto1in2 = (i < 3) ? i : ((i > 3) ? 2 : 1);
				break;
			}
		}
		for (i = 5; i >= 0; i--)
		{
			if (equalVerts[i])
			{
				eqto2in1 = (i < 3) ? 0 : ((i < 5) ? 1 : 2);
				eqto2in2 = (i < 3) ? i : ((i > 3) ? 2 : 1);
				break;
			}
		}
		p1 = verts1[eqto1in1];
		p2 = verts1[eqto2in1];
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
		p3 = verts1[lone1];
		p4 = verts2[lone2];

		// Compute cosines of angles p3-p2-p4 and p3-p1-p4
		Vec p3_p2 = Point.Difference (p3.Coord (), p2.Coord ());
		Vec p4_p2 = Point.Difference (p4.Coord (), p2.Coord ());
		Vec p3_p1 = Point.Difference (p3.Coord (), p1.Coord ());
		Vec p4_p1 = Point.Difference (p4.Coord (), p1.Coord ());
		double cos324 = p3_p2.Dot (p4_p2) /
		                (p3_p2.Length () * p4_p2.Length ());
		double cos314 = p3_p1.Dot (p4_p1) /
		                (p3_p1.Length () * p4_p1.Length ());

		// Return score = 1 - (1/2) |cos (p3-p2-p4) + cos (p3-p1-p4)|
		return 1.0 - 0.5 * Math.abs (cos324 + cos314);
	}

	// Sorts scores[] and correspondingly sorts pairs[]
	public void PairScoreQuicksort (double[] scores, TrianglePair[] pairs)
	{
		PairQuickSort (scores, pairs, 0, scores.length - 1);
	}

    void PairQuickSort (double[] qq, TrianglePair[] pp, int low, int high) {
		int lo = low;
		int hi = high;
		if (lo >= hi) {
			return;
		}
		double mid = qq [(lo+hi)/2];
		while (lo <= hi) {
			while (lo<high && mid < qq [lo]) {
				lo++;
			}
			while (hi>low && mid > qq [hi]) {
				hi--;
			}
			if (lo <= hi) {
				double w = qq [lo];
				qq [lo] = qq [hi];
				qq [hi] = w;
				TrianglePair pair = pp[lo];
				pp[lo] = pp[hi];
				pp[hi] = pair;
				lo++;
				hi--;
			}
		}
		
		if( low < hi )
            PairQuickSort (qq, pp, low, hi );
		if( lo < high )
            PairQuickSort (qq, pp, lo, high );
    }
};




