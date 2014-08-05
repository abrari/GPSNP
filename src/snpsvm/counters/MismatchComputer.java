package snpsvm.counters;

import java.util.Iterator;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FastaWindow;
import snpsvm.bamreading.FeatureComputer;
import snpsvm.bamreading.MappedRead;

public class MismatchComputer implements FeatureComputer {

    final double[] values = new double[1]; // Only compute ALT bases

    @Override
	public String getName(int which) {
		return "mismatch.alt";
	}
	
	@Override
	public int getColumnCount() {
		return values.length;
	}

	@Override
	public String getColumnDesc(int which) {
	    return "Number of mismatching bases on non-reference reads";
	}

	@Override
	public double[] computeValue(final char refBase, FastaWindow window, AlignmentColumn col) {
		values[0] = 0.0;

		if (col.getDepth() > 0) {
			Iterator<MappedRead> it = col.getIterator();
			while(it.hasNext()) {
				MappedRead read = it.next();
				if (read.hasBaseAtReferencePos(col.getCurrentPosition())) {
					char b = (char)read.getBaseAtReferencePos(col.getCurrentPosition());
					if (b == 'N' || refBase == 'N') 
						continue;
					
					int q = read.getMismatchCount(window);
                    // Only count variant base
					if ( b != refBase) {
						values[0] += q;
					}
				}
			}
		}

		return values;
	}

}
