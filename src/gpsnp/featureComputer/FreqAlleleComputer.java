package gpsnp.featureComputer;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FastaWindow;
import snpsvm.bamreading.MappedRead;

import java.util.Iterator;

/**
 * Frequency of major and minor allele
 *
 * Created by abrari on 4/28/14.
 */
public class FreqAlleleComputer extends AlleleComputer {
    @Override
    public String getName(int which) {
        if (which == major)
            return "freq.major";
        else
            return "freq.minor";
    }

    @Override
    public String getColumnDesc(int which) {
        if (which == major)
            return "Frequency of major allele";
        else
            return "Frequency of minor allele";
    }

    @Override
    public double[] computeValue(char refBase, FastaWindow window, AlignmentColumn col) {

        values[major] = 0.0;
        values[minor] = 0.0;

        int countMajor = 0;
        int countMinor = 0;

        if (col.getDepth() > 0) {

            calculateAlleles(col.getBases(), col.getDepth());
            Iterator<MappedRead> it = col.getIterator();

            while(it.hasNext()) {
                MappedRead read = it.next();
                if (read.hasBaseAtReferencePos(col.getCurrentPosition())) {
                    byte b = read.getBaseAtReferencePos(col.getCurrentPosition());
                    byte q = read.getQualityAtReferencePos(col.getCurrentPosition());

                    if ((char) b == getMajorAllele()) {
                        countMajor++;
                    } else if ((char) b == getMinorAllele()) {
                        countMinor++;
                    }
                }
            }

            values[major] = (double) countMajor / (double) col.getDepth();
            values[minor] = (double) countMinor / (double) col.getDepth();
        }

        return values;
    }
}
