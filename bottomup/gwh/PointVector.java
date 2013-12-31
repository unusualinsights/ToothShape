
import java.lang.*;

class Point
{
	double x, y, z;

	public Point (double x, double y, double z)
	{
		this.x = x;
		this.y = y;
		this.z = z;
	}

	public void Translate (Vec v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
	}

	public static Point Translate (Point p, Vec v)
	{
		return new Point (p.x + v.x, p.y + v.y, p.z + v.z);
	}

	public static Vec Difference (Point p, Point q)
	{
		return new Vec (p.x - q.x, p.y - q.y, p.z - q.z);
	}

	public double GetX () { return x; }
	public double GetY () { return y; }
	public double GetZ () { return z; }

	public String toString ()
	{
		StringBuffer sb = new StringBuffer ();
		sb.append ("(");
		sb.append (x);
		sb.append (", ");
		sb.append (y);
		sb.append (", ");
		sb.append (z);
		sb.append (")");
		return sb.toString ();
	}
};

class Vec
{
	double x, y, z;

	public Vec ()
	{
		x = 0.0; y = 0.0; z = 0.0;
	}

	public Vec (double x, double y, double z)
	{
		this.x = x;
		this.y = y;
		this.z = z;
	}

	public double Length ()
	{
		return Math.sqrt (x * x + y * y + z * z);
	}

	public void Add (Vec v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
	}

	public static Vec Add (Vec u, Vec v)
	{
		return new Vec (u.x + v.x, u.y + v.y, u.z + v.z);
	}

	public void Negate ()
	{
		x = -x;
		y = -y;
		z = -z;
	}

	public Vec Negation ()
	{
		return new Vec (-x, -y, -z);
	}

	public static Vec Negate (Vec v)
	{
		return new Vec (-v.x, -v.y, -v.z);
	}

	public void Normalize ()
	{
		double one_over_len = 1.0 / Length ();
		x *= one_over_len;
		y *= one_over_len;
		z *= one_over_len;
	}

	public static Vec Normalize (Vec v)
	{
		double one_over_len = 1.0 / v.Length ();
		return new Vec (v.x * one_over_len,
		                   v.y * one_over_len,
		                   v.z * one_over_len);
	}

	public void Subtract (Vec v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
	}

	public static Vec Subtract (Vec u, Vec v)
	{
		return new Vec (u.x - v.x, u.y - v.y, u.z - v.z);
	}

	public double Dot (Vec v)
	{
		return (x * v.x + y * v.y + z * v.z);
	}

	public static double Dot (Vec u, Vec v)
	{
		return (u.x * v.x + u.y * v.y + u.z * v.z);
	}

	public Vec Cross (Vec v)
	{
		return new Vec (y * v.z - z * v.y,
		                   z * v.x - x * v.z,
		                   x * v.y - y * v.x);
	}

	public static Vec Cross (Vec u, Vec v)
	{
		return new Vec (u.y * v.z - u.z * v.y,
		                   u.z * v.x - u.x * v.z,
		                   u.x * v.y - u.y * v.x);
	}

	public double GetX () { return x; }
	public double GetY () { return y; }
	public double GetZ () { return z; }

	public String toString ()
	{
		StringBuffer sb = new StringBuffer ();
		sb.append ("(");
		sb.append (x);
		sb.append (", ");
		sb.append (y);
		sb.append (", ");
		sb.append (z);
		sb.append (")");
		return sb.toString ();
	}
};

class Normal extends Vec
{
	public Normal (double x, double y, double z)
	{
		super (x, y, z);
	}
};


