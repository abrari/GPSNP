package gpsnp.featureComputer;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FastaWindow;
import snpsvm.bamreading.FeatureComputer;
import snpsvm.bamreading.MappedRead;

import java.util.Iterator;

/**
 * Mean base quality of ALL base in variant site
 * Not differentiating between ref/alt base
 *
 * Created by abrari on 4/28/14.
 */
public class BaseQualityComputer implements FeatureComputer {

    final double[] values = new double[1];

    @Override
    public String getName(int which) {
        return "mean.base.qual";
    }

    @Override
    public int getColumnCount() {
        return values.length;
    }

    @Override
    public String getColumnDesc(int which) {
        return "Mean base quality of reads in variant site";
    }

    @Override
    public double[] computeValue(char refBase, FastaWindow window, AlignmentColumn col) {
        values[0] = 0;
        int count = 0;

        if (col.getDepth() > 0) {
            Iterator<MappedRead> it = col.getIterator();
            while (it.hasNext()) {
                MappedRead read = it.next();
                if (read.hasBaseAtReferencePos(col.getCurrentPosition())) {
                    // Exclude ambiguous base
                    if (read.getBaseAtReferencePos(col.getCurrentPosition()) == 'N')
                        continue;

                    // System.err.print(read.getQualityAtReferencePos(col.getCurrentPosition()) + " ");

                    values[0] += read.getQualityAtReferencePos(col.getCurrentPosition());
                    count++;
                }
            }
        }

        // System.err.println();

        values[0] = values[0] / (double)count;

        return values;
    }
}
