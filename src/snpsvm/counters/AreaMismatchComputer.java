package snpsvm.counters;

import java.util.Iterator;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FeatureComputer;
import snpsvm.bamreading.FastaWindow;
import snpsvm.bamreading.MappedRead;

/**
 * Reports mean number of mismatches for all reads in window
 * @author brendan
 *
 */
public class AreaMismatchComputer implements FeatureComputer {

	private final double[] value = new double[1];
	
	@Override
	public String getName(int which) {
		return "area.mismatch";
	}
	
	@Override
	public double[] computeValue(final char refBase, FastaWindow window, AlignmentColumn col) {
		value[0] = 0;
        double counted = 0;
		if (col.getDepth() > 0) {
			Iterator<MappedRead> it = col.getIterator();
			while(it.hasNext()) {
				MappedRead read = it.next();
				if (read.hasBaseAtReferencePos(col.getCurrentPosition())) {
					char b = (char)read.getBaseAtReferencePos(col.getCurrentPosition());
					if (b == 'N' || refBase == 'N') 
						continue;
					
					value[0] += read.getMismatchCount(window);
                    counted++;
				}
			}
		}
		
		if (counted > 0)
			value[0] = value[0] / counted;

		return value;
	}

	@Override
	public int getColumnCount() {
		return value.length;
	}


	@Override
	public String getColumnDesc(int which) {
		return "Mean number of mismatches per read in window";
	}

}
