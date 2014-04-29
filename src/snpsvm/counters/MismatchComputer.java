package snpsvm.counters;

import java.util.Iterator;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FastaWindow;
import snpsvm.bamreading.MappedRead;

public class MismatchComputer extends VarCountComputer {
	
	@Override
	public String getName(int which) {
        if (which == ref)
		    return "mismatch.ref";
        else
            return "mismatch.alt";
	}
	
	@Override
	public int getColumnCount() {
		return values.length;
	}


	@Override
	public String getColumnDesc(int which) {
		if (which == ref)
			return "Number of mismatching bases on reference reads";
		else
			return "Number of mismatching bases on non-reference reads";
	}

	@Override
	public double[] computeValue(final char refBase, FastaWindow window, AlignmentColumn col) {
		values[ref] = 0.0;
		values[alt] = 0.0;

		if (col.getDepth() > 0) {
			Iterator<MappedRead> it = col.getIterator();
			while(it.hasNext()) {
				MappedRead read = it.next();
				if (read.hasBaseAtReferencePos(col.getCurrentPosition())) {
					char b = (char)read.getBaseAtReferencePos(col.getCurrentPosition());
					if (b == 'N' || refBase == 'N') 
						continue;
					
					int q = read.getMismatchCount(window);
					if ( b != refBase) {
						values[alt] += q;
					} else {
                        values[ref] += q;
                    }
				}
			}
		}

		return values;
	}

}
