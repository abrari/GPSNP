package snpsvm.counters;

import java.util.Iterator;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FeatureComputer;
import snpsvm.bamreading.FastaWindow;
import snpsvm.bamreading.MappedRead;

/*
 * Transition/transversion
 *
 * Transition = 0
 * Transversion = 1
 */

public class TsTvComputer implements FeatureComputer {

	final int[] counts = new int[4]; //Stores counts of each observed base
	final double[] values = new double[1]; //Return value
	
	@Override
	public String getName(int which) {
		return "ts.tv";
	}

	@Override
	public int getColumnCount() {
		return values.length;
	}

	@Override
	public String getColumnDesc(int which) {
		return "Mutation type, 0 if Ts, 1 if Tv";
	}

	@Override
	public double[] computeValue(char refBase, FastaWindow window,
			AlignmentColumn col) {

		counts[0] = 0;
		counts[1] = 0;
		counts[2] = 0;
		counts[3] = 0;
		
		int totAlts = 0;
		double val = 0.0;
		
		values[0] = 0;
		
		if (col.getDepth() > 0) {
			Iterator<MappedRead> it = col.getIterator();
			while(it.hasNext()) {
				MappedRead read = it.next();
				if (read.hasBaseAtReferencePos(col.getCurrentPosition())) {
					byte b = read.getBaseAtReferencePos(col.getCurrentPosition());
					
					int index = -1;
					if (b != refBase && b != 'N') {
						switch(b) {
						case 'A' : index = 0;
						case 'C' : index = 1;
						case 'G' : index = 2;
						case 'T' : index = 3;
						}
						counts[index]++;
						totAlts++;
					}
					
				}
			}
			
			char alt = computeAlt(counts);
			
			if ((refBase == 'A' && alt=='G') 
				|| (refBase == 'G' && alt=='A')
				|| (refBase == 'C' && alt=='T')
				|| (refBase == 'T' && alt=='C')) {
				val = 0; // Transition
			}
			else {
				val = 1; // Transversion
			}
			
			
		}
		
		values[0] = val;
		return values;
	}

	
	protected static int toIndex(char base) {
		if (base == 'A')
			return 0;
		if (base == 'C')
			return 1;
		if (base == 'G')
			return 2;
		if (base == 'T')
			return 3;
		return -1;
	}
	
	protected static char computeAlt(int[] counts) {
		int A = counts[0];
		int C = counts[1];
		int G = counts[2];
		int T = counts[3];
		
		//Find max of all...
		if (A >= C && A>=G && A>=T) {
			return 'A';
		}
		if (C >= A && C>=G && C>=T) {
			return 'C';
		}
		if (G >= A && G>=C && G>=T) {
			return 'G';
		}
		if (T >= A && T>=C && T>=G) {
			return 'T';
		}


		return '?';
	}
}
