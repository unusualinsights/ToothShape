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

	public Vector FuzzySegmentation (int k)
	{
		int N = trimesh.nfaces;
		Hashtable segments = new Hashtable ();

		double[] distances = new double[N];
		double minsumdist = Double.MAX_VALUE;
		int minindex = 0;
		for (int i = 0; i < N; i++)
		{
			GetDistances (i, distances);
if (i % 100 == 0)
System.err.println ("Distances for face " + i + " computed.");
			double sumdist = 0.0;
			for (int j = 0; j < N; j++)
				sumdist += distances[j];
			if (sumdist < minsumdist)
			{
				minsumdist = sumdist;
				minindex = i;
			}
		}

		// Indices of representative faces
		int[] rep = new int[k];
		rep[0] = minindex;
System.err.println ("rep[0] is face number " + rep[0]);

		// Faces not yet assigned as representatives
		Vector faces = new Vector ();
		for (int i = 0; i < N; i++)
			if (i != rep[0])
				faces.addElement (new Integer (i));

		// Compute indices of representative faces 1..k-1
		for (int j = 1; j < k; j++)
		{
			// Find face i w/ maximum min distance from other reps
			double maxmindist = 0.0;
			int maxminindex = 0;
			for (int f = 0; f < faces.size (); f++)
			{
				int i =
				    ((Integer)faces.elementAt (f)).intValue ();
				GetDistances (i, distances);
if (i % 100 == 0)
System.err.println ("j = " + j + ", distances for face " + i + " computed.");
				double mindist = distances[rep[0]];
				for (int m = 1; m < j; m++)
				{
					double imdist = distances[rep[m]];
					if (imdist < mindist)
						mindist = imdist;
				}
				if (mindist > maxmindist)
				{
					maxmindist = mindist;
					maxminindex = i;
				}
			}
			rep[j] = maxminindex;
System.err.println ("rep[" + j + " ] = " + rep[j]);
			faces.removeElement (new Integer (maxminindex));
		}

		// Compute probabilities[i][j], face i, patch j
		double[][] probs = new double[N][k];
		for (int i = 0; i < N; i++)
		{
			GetDistances (i, distances);
if (i % 100 == 0)
System.err.println ("Prob; distances for face " + i + " computed.");
			for (int j = 0; j < k; j++)
			{
				if (i == rep[j])
					probs[i][j] = 1.0;
				else
				{
					double numerator = 1.0 /
					                 distances[rep[j]];
					double denomsum = 0.0;
					for (int l = 0; l < k; l++)
						denomsum += 1.0 /
						         distances[rep[l]];
					probs[i][j] = numerator / denomsum;
				}
			}
		}

		// Nonfuzzy addition of each face to most likely segment
		for (int i = 0; i < N; i++)
		{
			Integer I = new Integer (i);
			double maxprob = 0.0;
			int maxsegnum = 0;
			for (int j = 0; j < k; j++)
			{
				if (probs[i][j] > maxprob)
				{
					maxprob = probs[i][j];
					maxsegnum = j;
				}
			}
			Integer Segnum = new Integer (maxsegnum);
			segments.put (I, Segnum);
		}

		Vector segInfo = new Vector ();
		segInfo.addElement (segments);
		segInfo.addElement (new Integer (k));
		segInfo.addElement (new Integer (N));

		return segInfo;
	}

	public Vector FuzzySegmentationForSmallSets (int k)
	{
		int N = trimesh.nfaces;
		Hashtable segments = new Hashtable ();

		// Matrix distances[i][j] is distance between faces i, j
		// Yes this is wasteful, distances[i][j] = distances[j][i]
		double[][] distances = new double[N][N];
		GetAllDistances (distances);

		double minsumdist = Double.MAX_VALUE;
		int minindex = 0;
		for (int i = 0; i < N; i++)
		{
			double sumdist = 0.0;
			for (int j = 0; j < N; j++)
				sumdist += distances[i][j];
			if (sumdist < minsumdist)
			{
				minsumdist = sumdist;
				minindex = i;
			}
		}

		// Indices of representative faces
		int[] rep = new int[k];
		rep[0] = minindex;

		// Faces not yet assigned as representatives
		Vector faces = new Vector ();
		for (int i = 0; i < N; i++)
			if (i != rep[0])
				faces.addElement (new Integer (i));

		// Compute indices of representative faces 1..k-1
		for (int j = 1; j < k; j++)
		{
			// Find face i w/ maximum min distance from other reps
			double maxmindist = 0.0;
			int maxminindex = 0;
			for (int f = 0; f < faces.size (); f++)
			{
				int i =
				    ((Integer)faces.elementAt (f)).intValue ();
				double mindist = distances[i][rep[0]];
				for (int m = 1; m < j; m++)
				{
					double imdist = distances[i][rep[m]];
					if (imdist < mindist)
						mindist = imdist;
				}
				if (mindist > maxmindist)
				{
					maxmindist = mindist;
					maxminindex = i;
				}
			}
			rep[j] = maxminindex;
			faces.removeElement (new Integer (maxminindex));
		}

		// Compute probabilities[i][j], face i, patch j
		double[][] probs = new double[N][k];
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < k; j++)
			{
				if (i == rep[j])
					probs[i][j] = 1.0;
				else
				{
					double numerator = 1.0 /
					                 distances[i][rep[j]];
					double denomsum = 0.0;
					for (int l = 0; l < k; l++)
						denomsum += 1.0 /
						         distances[i][rep[l]];
					probs[i][j] = numerator / denomsum;
				}
			}
		}

		// Nonfuzzy addition of each face to most likely segment
		for (int i = 0; i < N; i++)
		{
			Integer I = new Integer (i);
			double maxprob = 0.0;
			int maxsegnum = 0;
			for (int j = 0; j < k; j++)
			{
				if (probs[i][j] > maxprob)
				{
					maxprob = probs[i][j];
					maxsegnum = j;
				}
			}
			Integer Segnum = new Integer (maxsegnum);
			segments.put (I, Segnum);
		}

		Vector segInfo = new Vector ();
		segInfo.addElement (segments);
		segInfo.addElement (new Integer (k));
		segInfo.addElement (new Integer (N));

		return segInfo;
	}

	// Using Dijkstra's Algorithm, computes minimum distances between all
	// pairs of triangles using the dual graph of the triangle mesh, where
	// each vertex is the centroid of a triangle and each edge is the
	// distance between centroids of adjacent triangles
	public void GetAllDistances (double[][] distances)
	{
		int N = trimesh.nfaces;

		// Compute edgeshare array
		int[][] edgeshare = new int[N][3];
		for (int i = 0; i < N; i++)
		{
			int curIndex = 0;
			for (int j = 0; j < N; j++)
			{
				if (i != j && ShareEdge (i, j))
				{
					edgeshare[i][curIndex] = j;
					curIndex++;
				}
			}
		}

		for (int i = 0; i < N; i++)
			GetDijkstra (i, distances, edgeshare);
	}

	// Compute shortest distances from faces[index] to all other faces
	public void GetDijkstra (int index, double[][] dist, int[][] neighbors)
	{
		int N = trimesh.nfaces;
		int i;

		for (i = 0; i < N; i++)
			dist[index][i] = (i == index) ? 0 : Double.MAX_VALUE;

		Vector S = new Vector ();
		Vector Q = new Vector ();
		for (i = 0; i < N; i++)
			Q.addElement (new Integer (i));

		while (!Q.isEmpty ())
		{
			// Extract element in Q with minimum dist[] value
			int u = ((Integer)Q.elementAt (0)).intValue ();
			for (i = 1; i < Q.size (); i++)
			{
				int q = ((Integer)Q.elementAt (i)).intValue ();
				if (dist[index][q] < dist[index][u]) u = q;
			}
			Integer U = new Integer (u);
			Q.removeElement (U);

			S.addElement (U);
			// Relax edges incident on u
			for (i = 0; i < 3; i++)
			{
				int nbr = neighbors[u][i];
				// Relax edge (u, nbr)
				double dui = GetDistance (u, nbr);
				if (dist[index][nbr] > dist[index][u] + dui)
					dist[index][nbr] = dist[index][u] + dui;
			}
		}
	}

	// Returns Dijkstra-based distance from triangle index to all other
	// triangles
	public void GetDistances (int index, double[] dist)
	{
		int N = trimesh.nfaces;
		GetDijkstra (index, dist);
	}

	public void GetDijkstra (int index, double[] dist)
	{
		int N = trimesh.nfaces;
		int i;

		for (i = 0; i < N; i++)
			dist[i] = (i == index) ? 0 : Double.MAX_VALUE;

		Vector S = new Vector ();
		Vector Q = new Vector ();
		for (i = 0; i < N; i++)
			Q.addElement (new Integer (i));

		while (!Q.isEmpty ())
		{
			// Extract element in Q with minimum dist[] value
			int u = ((Integer)Q.elementAt (0)).intValue ();
			for (i = 1; i < Q.size (); i++)
			{
				int q = ((Integer)Q.elementAt (i)).intValue ();
				if (dist[q] < dist[u]) u = q;
			}
			Integer U = new Integer (u);
			Q.removeElement (U);

			S.addElement (U);
			// Relax edges incident on u
			for (i = 0; i < 3; i++)
			{
				int nbr = neighbors[u][i];
				// Relax edge (u, nbr)
				double dui = GetDistance (u, nbr);
				if (dist[nbr] > dist[u] + dui)
					dist[nbr] = dist[u] + dui;
			}
		}
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




