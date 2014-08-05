package snpsvm.counters;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FeatureComputer;
import snpsvm.bamreading.FastaWindow;

public class DinucRepeatCounter implements FeatureComputer {

	final double[] values = new double[1]; // now we combine left and right value
	
	final int maxLength = 10; //dont look beyond this many bases in either direction
	
	@Override
	public String getName(int which) {
        return "dinuc.repeat";
	}
	
	@Override
	public int getColumnCount() {
		return values.length;
	}


	@Override
	public String getColumnDesc(int which) {
		return "Number of Dinucleotide repeats to left and right of site";
	}

	@Override
	public double[] computeValue(char refBase, FastaWindow window, AlignmentColumn col) {

        char base0, base1;
        double count;
        double left, right;
		
		//Looking backward
		int refPos = col.getCurrentPosition();
        // Fix if near left edge of contig
        if(refPos <= 2) {
            left = 0.0;
        } else {
            base0 = window.getBaseAt(refPos-1);
            base1 = window.getBaseAt(refPos-2);
            count = 0;
            if (base0 != base1) {
                for(int i=refPos-3; (i-1)>Math.max(window.indexOfLeftEdge(), refPos-1-maxLength); i-=2) {
                    if (base0 == window.getBaseAt(i) && base1 == window.getBaseAt(i-1))
                        count++;
                    else
                        break;
                }
            }
            left = count;
        }
		
		//Looking forward
        // Fix if near right edge of window where refPos+2 will fail
        if(refPos >= window.indexOfRightEdge() - 2) {
            right = 0.0;
        } else {
            base0 = window.getBaseAt(refPos + 1);
            base1 = window.getBaseAt(refPos + 2);
            count = 0;
            if (base0 != base1) {
                for (int i = refPos + 3; (i + 1) < Math.min(window.indexOfRightEdge(), refPos + 1 + maxLength); i += 2) {
                    if (base0 == window.getBaseAt(i) && base1 == window.getBaseAt(i + 1))
                        count++;
                    else
                        break;
                }
            }
            right = count;
        }

		values[0] = left + right;
		return values;
	}
}
