package snpsvm.counters;

import java.util.Iterator;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FeatureComputer;
import snpsvm.bamreading.FastaWindow;
import snpsvm.bamreading.MappedRead;

public class NearbyQualComputer implements FeatureComputer {

    private boolean isPhred33;

	public final int WINDOW_SIZE = 5; //Window spans the focus position, so 7 means three in either direction
	double[] values = new double[WINDOW_SIZE];
	double[] counts = new double[WINDOW_SIZE];

    public NearbyQualComputer() {

    }

    public NearbyQualComputer(boolean isPhred33) {
        this.isPhred33 = isPhred33;
    }

	@Override
	public String getName(int which) {
		return "mean.nearby.qual";
	}

	@Override
	public int getColumnCount() {
		return 1;
	}


	@Override
	public String getColumnDesc(int which) {
		return "Mean quality of sites aligning nearby";
	}
	
	@Override
	public double[] computeValue(final char refBase, FastaWindow window, AlignmentColumn col) {
		for(int i=0; i<WINDOW_SIZE; i++) {
			values[i] = 0.0;
			counts[i] = 0.0;
		}
		
		int offset = WINDOW_SIZE/2;
		if (col.getDepth() > 0) {
			Iterator<MappedRead> it = col.getIterator();
			while(it.hasNext()) {
				MappedRead read = it.next();
				for(int i=0; i<WINDOW_SIZE; i++) {
					int refPos = col.getCurrentPosition()-offset+i;
                    if(refPos != col.getCurrentPosition()) { // Exclude current column
                        if (read.hasBaseAtReferencePos(refPos)) {
                            values[i] += read.getQualityAtReferencePos(refPos);
                            if (isPhred33) values[i] += 31;

                            counts[i]++;
                        }
                    }
                }
            }
		}

        double sumOfMeans = 0.0;
        double allQualMeans;

		for (int i=0; i<WINDOW_SIZE; i++) {
            if (i != WINDOW_SIZE/2) {
                if (counts[i] > 0) {
                    values[i] /= counts[i];
                }
                sumOfMeans += values[i];
            }
		}

        if (WINDOW_SIZE > 1)
            allQualMeans = sumOfMeans / (double)(WINDOW_SIZE-1);

		return new double[]{allQualMeans};
	}

}
