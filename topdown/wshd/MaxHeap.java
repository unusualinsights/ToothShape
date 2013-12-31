
class FMWData
{
	public int vertex;
	public boolean label;

	public FMWData (int v, boolean lbl)
	{
		vertex = v;
		label = lbl;
	}
};

class MaxHeap
{
	double[] keys;
	FMWData[] values;
	int numElements;

	public MaxHeap ()
	{
		numElements = 0;
	}

	public MaxHeap (double[] ks, FMWData[] vs)
	{
		numElements = ks.length;
		keys = new double[numElements];
		values = new FMWData[numElements];
		for (int i = 0; i < numElements; i++)
		{
			keys[i] = ks[i];
			values[i] = vs[i];
		}
		for (int i = numElements / 2 - 1; i >= 0; i--)
			MaxHeapify (i);
	}

	public int size ()
	{
		return numElements;
	}

	public boolean isEmpty ()
	{
		return numElements == 0;
	}

	public FMWData Pop ()
	{
		if (numElements < 1)
		{
			System.err.println ("Heap is empty.");
			System.exit (0);
		}
		double max = keys[0]; // this value is useless, not returned
		FMWData maxval = values[0]; // this value is returned, can get
                                            // max which is k2[maxval.vertex]
		keys[0] = keys[numElements - 1];
		values[0] = values[numElements - 1];
		numElements--;
		MaxHeapify (0);
		return maxval;
	}

	public void Push (double key, FMWData value)
	{
		if (numElements > 0)
		{
			double[] temp = new double[numElements];
			FMWData[] fmwTmp = new FMWData[numElements];
			for (int i = 0; i < numElements; i++)
			{
				temp[i] = keys[i];
				fmwTmp[i] = values[i];
			}
			keys = new double[numElements + 1];
			values = new FMWData[numElements + 1];
			for (int i = 0; i < numElements; i++)
			{
				keys[i] = temp[i];
				values[i] = fmwTmp[i];
			}
			keys[numElements] = Double.NEGATIVE_INFINITY;
			numElements++;
			IncreaseKey (numElements - 1, key, value);
		}
		else
		{
			keys = new double[1];
			values = new FMWData[1];
			keys[0] = Double.NEGATIVE_INFINITY;
			numElements = 1;
			IncreaseKey (0, key, value);
		}
	}

	public void IncreaseKey (int i, double key, FMWData value)
	{
		if (key < keys[i])
		{
			System.err.println ("New key smaller than current");
			System.exit (0);
		}
		keys[i] = key;
		values[i] = value;
		while (i > 0 && keys[Parent (i)] < keys[i])
		{
			double temp = keys[i];
			keys[i] = keys[Parent (i)];
			keys[Parent (i)] = temp;

			FMWData fmwTmp = values[i];
			values[i] = values[Parent (i)];
			values[Parent (i)] = fmwTmp;

			i = Parent (i);
		}
	}

	public void MaxHeapify (int i)
	{
		if (i < 0 || i > numElements - 1)
			return;

		int l = Left (i);
		int r = Right (i);
		int largest = -1;
		if (l < numElements && keys[l] > keys[i])
			largest = l;
		else
			largest = i;
		if (r < numElements && keys[r] > keys[largest])
			largest = r;
		if (largest != i)
		{
			double temp = keys[i];
			keys[i] = keys[largest];
			keys[largest] = temp;

			FMWData fmwTmp = values[i];
			values[i] = values[largest];
			values[largest] = fmwTmp;

			MaxHeapify (largest);
		}
	}

	public double GetParent (int i)
	{
		return keys[Parent (i)];
	}

	public double GetLeft (int i)
	{
		return keys[Left (i)];
	}

	public double GetRight (int i)
	{
		return keys[Right (i)];
	}

	public static int Parent (int i)
	{
		return (i - 1) / 2;
	}

	public static int Left (int i)
	{
		return 2 * i + 1;
	}

	public static int Right (int i)
	{
		return Left (i) + 1;
	}
};


