import MeshSegmenter;
import java.lang.*;
import java.util.*;
import java.util.Vector;

public class pairtest
{
	public static void main (String args[])
	{
		//TriangleMesh trimesh = new TriangleMesh ("tstmesh.jtm");
		TriangleMesh trimesh = new TriangleMesh ("../../jdata/CaseCFPREcut_Tooth_19.jtm");
		MeshSegmenter ms = new MeshSegmenter (trimesh);
		TrianglePair[] pairs = ms.GetPairs ();
		int i = 0;
		Vertex[] verts = ms.Vertices (pairs[i]);
		Face tri1 = trimesh.faces[pairs[i].face1];
		Face tri2 = trimesh.faces[pairs[i].face2];
		Vertex p10 = trimesh.vertices[tri1.fvlist[0]];
		Vertex p11 = trimesh.vertices[tri1.fvlist[1]];
		Vertex p12 = trimesh.vertices[tri1.fvlist[2]];
		Vertex p20 = trimesh.vertices[tri2.fvlist[0]];
		Vertex p21 = trimesh.vertices[tri2.fvlist[1]];
		Vertex p22 = trimesh.vertices[tri2.fvlist[2]];
		System.out.print ("Pair numbers: (");
		System.out.print (pairs[i].face1);
		System.out.print (", ");
		System.out.print (pairs[i].face2);
		System.out.println (")");
		System.out.println ("Vertices of face 1:");
		System.out.println (p10.toString ());
		System.out.println (p11.toString ());
		System.out.println (p12.toString ());
		System.out.println ("Vertices of face 2:");
		System.out.println (p20.toString ());
		System.out.println (p21.toString ());
		System.out.println (p22.toString ());
		System.out.println ("Same vertices:");
		System.out.println (verts[0].toString ());
		System.out.println (verts[1].toString ());
		System.out.println ("Different vertices:");
		System.out.println (verts[2].toString ());
		System.out.println (verts[3].toString ());
/*
		for (int i = 0; i < pairs.length; i++)
		{
			StringBuffer sb = new StringBuffer ();
			sb.append ("(");
			sb.append (pairs[i].face1);
			sb.append (", ");
			sb.append (pairs[i].face2);
			sb.append (")");
			System.out.println (sb.toString ());
		}
		Vec norm1 = ms.UnnormalizedFaceNormal (trimesh.faces[2]);
		Vec norm2 = ms.UnnormalizedFaceNormal (trimesh.faces[3]);
		Vec norm3 = ms.UnnormalizedFaceNormal (trimesh.faces[0]);
		norm1.Normalize ();
		norm2.Normalize ();
		norm3.Normalize ();
		System.out.print ("Face 2 normal: (");
		System.out.print (norm1.GetX () + ", ");
		System.out.print (norm1.GetY () + ", ");
		System.out.println (norm1.GetZ () + ")");
		System.out.print ("Face 3 normal: (");
		System.out.print (norm2.GetX () + ", ");
		System.out.print (norm2.GetY () + ", ");
		System.out.println (norm2.GetZ () + ")");
		System.out.print ("Face 0 normal: (");
		System.out.print (norm3.GetX () + ", ");
		System.out.print (norm3.GetY () + ", ");
		System.out.println (norm3.GetZ () + ")");
		System.out.println ("2 dot 3 = " + norm1.Dot (norm2));
		System.out.println ("2 dot 0 = " + norm1.Dot (norm3));
*/
	}
};


