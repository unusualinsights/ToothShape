
import java.lang.*;
import java.io.*;
import java.util.*;
//import PointVector;

class Vertex
{
	Point coordinates;
	boolean has_normals;
	Normal normal;

	public Vertex (double x, double y, double z)
	{
		coordinates = new Point (x, y, z);
		has_normals = false;
	}

	public Vertex (double  x, double  y, double  z,
	               double nx, double ny, double nz)
	{
		coordinates = new Point (x, y, z);
		has_normals = true;
		normal = new Normal (nx, ny, nz);
	}

	public Point Coord ()
	{
		return coordinates;
	}

	public double GetX () { return coordinates.GetX (); }
	public double GetY () { return coordinates.GetY (); }
	public double GetZ () { return coordinates.GetZ (); }

	public String toString ()
	{
		StringBuffer sb = new StringBuffer ();
		sb.append ("(");
		sb.append (coordinates.GetX ());
		sb.append (", ");
		sb.append (coordinates.GetY ());
		sb.append (", ");
		sb.append (coordinates.GetZ ());
		sb.append (")");
		return sb.toString ();
	}
};

class Face
{
	public int[] fvlist;

	public Face ()
	{
		fvlist = new int[3];
	}
};

class TriangleMesh
{
	public boolean has_normals;
	public int nverts, nfaces;
	public Vertex[] vertices;
	public Face[] faces;

	public TriangleMesh () {}

	public TriangleMesh (String inJtmFileName)
	{
		try
		{
		BufferedReader in = new BufferedReader (new FileReader
		                                            (inJtmFileName));
		String line;
		int i;

		line = in.readLine ().trim ();
		if (line == null)
		{
			System.out.println ("Bad JTM file: " + inJtmFileName);
			System.out.println ("Was expecting has_normals");
			System.out.println ("Found end of file instead");
			System.exit (0);
		}
		has_normals = line.equals ("yes_have_normals");

		line = in.readLine ().trim ();
		if (line == null)
		{
			System.out.println ("Bad JTM file: " + inJtmFileName);
			System.out.println ("Was expecting number of vertices");
			System.out.println ("Found end of file instead");
			System.exit (0);
		}
		nverts = Integer.parseInt (line);
		vertices = new Vertex[nverts];

		line = in.readLine ().trim ();
		if (line == null)
		{
			System.out.println ("Bad JTM file: " + inJtmFileName);
			System.out.println ("Was expecting number of triangles");
			System.out.println ("Found end of file instead");
			System.exit (0);
		}
		nfaces = Integer.parseInt (line);
		faces = new Face[nfaces];

		for (i = 0; i < nverts; i++)
		{
			line = in.readLine ().trim ();
			if (line == null)
			{
				System.out.println ("Bad JTM file: " + inJtmFileName);
				System.out.println ("Expected number of vertices: " + nverts);
				System.out.println ("Found only " + i + " vertices.");
				System.out.println ("And file ended unexpectedly.");
				System.exit (0);
			}

			StringTokenizer st = new StringTokenizer (line);
			double px = Double.parseDouble (st.nextToken ());
			double py = Double.parseDouble (st.nextToken ());
			double pz = Double.parseDouble (st.nextToken ());
			if (has_normals)
			{
				double nx = Double.parseDouble (st.nextToken ());
				double ny = Double.parseDouble (st.nextToken ());
				double nz = Double.parseDouble (st.nextToken ());
				vertices[i] = new Vertex (px, py, pz, nx, ny, nz);
			}
			else
				vertices[i] = new Vertex (px, py, pz);
		}

		for (i = 0; i < nfaces; i++)
		{
			line = in.readLine ().trim ();
			if (line == null)
			{
				System.out.println ("Bad JTM file: " + inJtmFileName);
				System.out.println ("Expected number of triangles: " + nfaces);
				System.out.println ("Found only " + i + " faces.");
				System.out.println ("And file ended unexpectedly.");
				System.exit (0);
			}

			StringTokenizer st = new StringTokenizer (line);
			faces[i] = new Face ();
			faces[i].fvlist[0] = Integer.parseInt (st.nextToken ());
			faces[i].fvlist[1] = Integer.parseInt (st.nextToken ());
			faces[i].fvlist[2] = Integer.parseInt (st.nextToken ());
		}

		in.close ();
		}
		catch (Exception e) {e.printStackTrace ();}
	}

	void OutputJTM ()
	{
		int i;

		if (has_normals)
			System.out.println ("yes_have_normals");
		else
			System.out.println ("not_have_normals");

		System.out.println (nverts);
		System.out.println (nfaces);

		for (i = 0; i < nverts; i++)
		{
			System.out.print (vertices[i].coordinates.x + " ");
			System.out.print (vertices[i].coordinates.y + " ");
			System.out.print (vertices[i].coordinates.z + " ");
			System.out.print (vertices[i].normal.x + " ");
			System.out.print (vertices[i].normal.y + " ");
			System.out.print (vertices[i].normal.z);
			System.out.println ();
		}

		for (i = 0; i < nfaces; i++)
		{
			System.out.print (faces[i].fvlist[0] + " ");
			System.out.print (faces[i].fvlist[1] + " ");
			System.out.print (faces[i].fvlist[2]);
			System.out.println ();
		}
	}

	public TriangleMesh FaceToTriangleMesh (int i)
	{
		TriangleMesh face = new TriangleMesh ();

		face.has_normals = has_normals;
		face.nverts = 3;
		face.nfaces = 1;
		face.vertices = new Vertex[3];
		face.vertices[0] = vertices[faces[i].fvlist[0]];
		face.vertices[1] = vertices[faces[i].fvlist[1]];
		face.vertices[2] = vertices[faces[i].fvlist[2]];
		face.faces = new Face[1];
		face.faces[0].fvlist[0] = faces[i].fvlist[0];
		face.faces[0].fvlist[1] = faces[i].fvlist[1];
		face.faces[0].fvlist[2] = faces[i].fvlist[2];

		return face;
	}
};




