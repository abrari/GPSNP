package snpsvm.counters;

import java.util.Iterator;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FeatureComputer;
import snpsvm.bamreading.FastaWindow;
import snpsvm.bamreading.MappedRead;
import snpsvm.util.BinomMath;

public class ErrorProbComputer implements FeatureComputer {

	double[] value = new double[1];
	
	@Override
	public int getColumnCount() {
		return value.length;
	}


	@Override
	public String getColumnDesc(int which) {
		return "Probability base was not sampled from binomial distribution with p=0.05";
	}
	
	@Override
	public double[] computeValue(final char refBase, FastaWindow window, AlignmentColumn col) {
		int refCount = 0;
		int altCount = 0;
		if (col.getDepth() > 0) {
			Iterator<MappedRead> it = col.getIterator();
			while(it.hasNext()) {
				MappedRead read = it.next();
				if (read.hasBaseAtReferencePos(col.getCurrentPosition())) {
					byte b = read.getBaseAtReferencePos(col.getCurrentPosition());
					if (b == 'N') 
						continue;
					if (b == refBase)
						refCount++;
					else
						altCount++;
				}
			}
		}
		
		
		
		int T = refCount + altCount;
		int X = altCount;

        // Cap at 250
		if (T > 250) {
			X = (250 * X)/T;
			T = 250;
		}
		
		
		//Compute het prob
		//Each read has 50% chance of coming from source with a non-reference base
		double hetProb = BinomMath.binomPDF((int)Math.round(X), (int) Math.round(T), 0.5);
		

		//Compute homo prob
		double homProb = BinomMath.binomPDF((int) Math.round(X), (int) Math.round(T), 0.99);
		
		//Compute error prob
        //Example:
        //   total=10, alt=1, then errorProb = 0.3
        //   total=10, alt=2, then errorProb = 0.07
		double errProb = BinomMath.binomPDF((int) Math.round(X), (int) Math.round(T), 0.05);

        value[0] = errProb / (hetProb + homProb + errProb);
		
		return value;
	}
	
	@Override
	public String getName(int which) {
		return "error.prob";
	}
	
}
