package gpsnp.featureComputer;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FastaWindow;
import snpsvm.bamreading.FeatureComputer;
import snpsvm.bamreading.MappedRead;

import java.util.Iterator;

/**
 * Mean mapping quality of ALL reads in variant site
 * Not differentiating between ref/alt base
 *
 * Created by abrari on 4/28/14.
 */
public class AlignmentQualityComputer implements FeatureComputer {

    final double[] values = new double[1];

    @Override
    public String getName() {
        return "mapping.qual";
    }

    @Override
    public int getColumnCount() {
        return values.length;
    }

    @Override
    public String getColumnDesc(int which) {
        return "Mean mapping quality of reads in variant site";
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

                    // System.err.print(read.getRecord().getMappingQuality() + " ");

                    values[0] += read.getRecord().getMappingQuality();
                    count++;
                }
            }
        }

        // System.err.println();

        values[0] = values[0] / (double)count;

        return values;
    }
}
