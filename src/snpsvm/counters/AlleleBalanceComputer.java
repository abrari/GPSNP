package snpsvm.counters;

import java.util.Iterator;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FeatureComputer;
import snpsvm.bamreading.FastaWindow;
import snpsvm.bamreading.MappedRead;

/**
 * Just emits fraction of reads with variant base
 * @author brendan
 *
 */
public class AlleleBalanceComputer implements FeatureComputer {
	
	double[] value = new double[1];
	
	@Override
	public double[] computeValue(final char refBase, FastaWindow window, AlignmentColumn col) {
		double count = 0;
		double varCount = 0;
		if (col.getDepth() > 0) {
			Iterator<MappedRead> it = col.getIterator();
			while(it.hasNext()) {
				MappedRead read = it.next();
				if (read.hasBaseAtReferencePos(col.getCurrentPosition())) {
					byte b = read.getBaseAtReferencePos(col.getCurrentPosition());
					if (b == 'N')
						continue;

					if (b != refBase)
						varCount++;

					count++;
				}
			}
		}
		
	    if (count > 0)
            value[0] = (double)varCount / (double)count;

		return value;
	}
	

	@Override
	public String getName(int which) {
		return "allele.balance";
	}


	@Override
	public int getColumnCount() {
		return value.length;
	}


	@Override
	public String getColumnDesc(int which) {
		return "Fraction of variant bases at site";
	}

}
