//import MeshSegmenter;
import java.lang.*;
import java.util.*;
import java.util.Vector;

public class kkm
{
	public static void main (String args[])
	{
		TriangleMesh trimesh = new TriangleMesh (args[0]);
		MeshSegmenter ms = new MeshSegmenter (trimesh);
		Vector segInfo = new Vector ();
		try
		{
			int k = Integer.parseInt (args[1]);
			if (k < 2)
			{
				System.err.println ("k set to be 2.");
				k = 2;
			}
			segInfo = ms.KkmSegmentation (k);
		}
		catch (Exception ee) {ee.printStackTrace ();}
		Hashtable segments = (Hashtable)segInfo.elementAt (0);
		Integer NumSegments = (Integer)segInfo.elementAt (1);
		int numSegments = NumSegments.intValue ();
		Hashtable segToTris = new Hashtable ();
		for (Enumeration e = segments.keys (); e.hasMoreElements ();)
		{
			Integer triNum = (Integer)e.nextElement ();
			Integer segNum = (Integer)segments.get (triNum);
			Vector tris = (Vector)segToTris.get (segNum);
			if (tris == null) tris = new Vector ();
			tris.addElement (triNum);
			segToTris.remove (segNum);
			segToTris.put (segNum, tris);
		}

		System.out.println (numSegments);
		System.out.println ();
		for (Enumeration e = segToTris.keys (); e.hasMoreElements ();)
		{
			Integer segNum = (Integer)e.nextElement ();
			Vector tris = (Vector)segToTris.get (segNum);
			int numtris = tris.size ();
			System.out.println (numtris);
			for (int i = 0; i < numtris; i++)
			{
				Integer Tnum = (Integer)tris.elementAt (i);
				System.out.println (Tnum.intValue ());
			}
			System.out.println ();
		}
	}
};


