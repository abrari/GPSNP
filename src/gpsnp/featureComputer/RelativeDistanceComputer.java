package gpsnp.featureComputer;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FastaWindow;
import snpsvm.bamreading.FeatureComputer;
import snpsvm.bamreading.MappedRead;

import java.util.Iterator;

/**
 * Relative distance to both ends of read
 *
 * Created by abrari on 4/29/14.
 */
public class RelativeDistanceComputer implements FeatureComputer {

    final double[] values = new double[1];

    @Override
    public String getName(int which) {
        return "rel.dist";
    }

    @Override
    public int getColumnCount() {
        return values.length;
    }

    @Override
    public String getColumnDesc(int which) {
        return "Relative distance of variant position in reads to both ends.";
    }

    @Override
    public double[] computeValue(char refBase, FastaWindow window, AlignmentColumn col) {

        values[0] = 0.0;
        int posInRead, readLength;
        int count = 0;
        double rd = 0;

        if (col.getDepth() > 0) {
            Iterator<MappedRead> it = col.getIterator();
            while (it.hasNext()) {
                MappedRead read = it.next();
                if (read.hasBaseAtReferencePos(col.getCurrentPosition())) {

                    byte b = read.getBaseAtReferencePos(col.getCurrentPosition());
                    posInRead = read.refPosToReadPos(col.getCurrentPosition());
                    readLength = read.getRecord().getReadLength();

                    // Only count variant base
                    if (b != refBase) {
                        double posRatio = (double)posInRead / (double)readLength;
                        if (posRatio <= 0.5)
                            rd += posRatio;
                        else
                            rd += (1 - posRatio);

                        count++;
                    }
                }
            }

            if (count > 0) {
                values[0] = rd / (double) count;
            }
        }

        return values;
    }
}
