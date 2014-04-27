package snpsvm.counters;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.ColumnComputer;
import snpsvm.bamreading.FastaWindow;

/**
 * Computes longest homopolymer run in both directions 
 * @author brendan
 *
 */
public class HomopolymerRunCounter implements ColumnComputer {

	final double[] values = new double[2];
	
	final int maxLength = 10; //dont look beyond this many bases in either direction
	
	@Override
	public String getName() {
		return "hrun.counter";
	}

	@Override
	public int getColumnCount() {
		return values.length;
	}


	@Override
	public String getColumnDesc(int which) {
		if (which == 0)
			return "Length of homopolymer run to left of site";
		else
			return "Length of homopolymer run to right of site";
	}
	
	@Override
	public double[] computeValue(char refBase, FastaWindow window,
			AlignmentColumn col) {

        char base;
        double count;
		
		//Looking backward
		int refPos = col.getCurrentPosition();

        // Fix if near left edge of contig where refPos-1 will fail
        if(refPos <= 1) {
            values[0] = 0.0;
        } else {
            base = window.getBaseAt(refPos - 1);
            count = 0;
            for (int i = refPos - 2; i > Math.max(window.indexOfLeftEdge(), refPos - 1 - maxLength); i--) {
                if (base == window.getBaseAt(i))
                    count++;
                else
                    break;
            }
            values[0] = count;
        }

        //Looking forward

        // Fix if near right edge of window where refPos+1 will fail
        if(refPos >= window.indexOfRightEdge() - 1) {
            values[1] = 0.0;
        } else {
            base = window.getBaseAt(refPos + 1);
            count = 0;
            for (int i = refPos + 2; i < Math.min(window.indexOfRightEdge() - 1, refPos + 1 + maxLength); i++) {
                if (i < window.indexOfRightEdge()) {
                    //System.out.println("requested pos : " + i + " right edge : " + window.indexOfRightEdge());
                    if (base == window.getBaseAt(i))
                        count++;
                    else
                        break;
                } else
                    break;
            }
            values[1] = count;
        }
		
		values[0] = values[0] / maxLength * 2.0 -1.0;
		values[1] = values[1] / maxLength * 2.0 -1.0;
		
		return values;
	}

}
