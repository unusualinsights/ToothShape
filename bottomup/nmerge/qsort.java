
/**
Extract words in a file system and organize that as a sorted list
*/


public class qsort {

	static int array [] = { 9, 1, 8, 3, 5, 6, 2, 0, 1, 1};

	public qsort () {
		int n = array.length;
		System.out.println ("Before sorting");
		for (int i=0; i<n; i++) {
			System.out.print (" "+array [i]);
		}
		System.out.println ();
		quickSort (array, 0, n-1);
		System.out.println ("After sorting");
		for (int i=0; i<n; i++) {
			System.out.print (" "+array [i]);
		}
		System.out.println ();
	}
					
    void quickSort (int qq [], int low, int high) {
		int lo = low;
		int hi = high;
		if (lo >= hi) {
			return;
		}
		int mid = qq [(lo+hi)/2];
		while (lo <= hi) {
			while (lo<high && mid < qq [lo]) {
				lo++;
			}
			while (hi>low && mid > qq [hi]) {
				hi--;
			}
			if (lo <= hi) {
				int w = qq [lo];
				qq [lo] = qq [hi];
				qq [hi] = w;
				lo++;
				hi--;
			}
		}
		
		if( low < hi )
            quickSort (qq, low, hi );
		if( lo < high )
            quickSort (qq, lo, high );
    }

	public static void main (String args []) {
		qsort q = new qsort ();
	}
};


