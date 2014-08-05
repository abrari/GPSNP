package snpsvm.counters;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FeatureComputer;
import snpsvm.bamreading.FastaWindow;

/**
 * Computes longest homopolymer run in both directions 
 * @author brendan
 *
 */
public class HomopolymerRunComputer implements FeatureComputer {

	final double[] values = new double[1]; // now we combine left and right value
	
	final int maxLength = 10; //dont look beyond this many bases in either direction
	
	@Override
	public String getName(int which) {
        return "homopolymer.length";
	}

	@Override
	public int getColumnCount() {
		return values.length;
	}


	@Override
	public String getColumnDesc(int which) {
		return "Length of homopolymer run on reference to left and right of site";
	}
	
	@Override
	public double[] computeValue(char refBase, FastaWindow window, AlignmentColumn col) {

        char base;
        double count;
        double left, right;
		
		//Looking backward
		int refPos = col.getCurrentPosition();

        // Fix if near left edge of contig where refPos-1 will fail
        if(refPos <= 1) {
            left = 0.0;
        } else {
            base = window.getBaseAt(refPos - 1);
            count = 0;
            for (int i = refPos - 2; i > Math.max(window.indexOfLeftEdge(), refPos - 1 - maxLength); i--) {
                if (base == window.getBaseAt(i))
                    count++;
                else
                    break;
            }
            left = count;
        }

        //Looking forward

        // Fix if near right edge of window where refPos+1 will fail
        if(refPos >= window.indexOfRightEdge() - 1) {
            right = 0.0;
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
            right = count;
        }

        values[0] = left + right;
		return values;
	}

}
