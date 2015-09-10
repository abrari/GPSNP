package gpsnp.featureComputer;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FastaWindow;
import snpsvm.bamreading.MappedRead;

import java.util.Iterator;

/**
 * Mean quality of major and minor allele
 *
 * Created by abrari on 4/28/14.
 */
public class MeanQualAlleleComputer extends AlleleComputer {

    public MeanQualAlleleComputer() {

    }

    public MeanQualAlleleComputer(boolean isPhred33) {
        this.isPhred33 = isPhred33;
    }

    @Override
    public String getName(int which) {
        if (which == major)
            return "mean.qual.major";
        else
            return "mean.qual.minor";
    }

    @Override
    public String getColumnDesc(int which) {
        if (which == major)
            return "Mean quality of major allele";
        else
            return "Mean quality of minor allele";
    }

    @Override
    public double[] computeValue(char refBase, FastaWindow window, AlignmentColumn col) {

        values[major] = 0.0;
        values[minor] = 0.0;

        int countMajor = 0;
        int countMinor = 0;

        if (col.getDepth() > 0) {

            calculateAlleles(col.getBases(), col.getDepth(), refBase);
            Iterator<MappedRead> it = col.getIterator();

            while(it.hasNext()) {
                MappedRead read = it.next();
                if (read.hasBaseAtReferencePos(col.getCurrentPosition())) {
                    byte b = read.getBaseAtReferencePos(col.getCurrentPosition());
                    byte q = read.getQualityAtReferencePos(col.getCurrentPosition());

                    if (isPhred33) q += 31;

                    if ((char) b == getMajorAllele()) {
                        values[major] += q;
                        countMajor++;
                    } else if ((char) b == getMinorAllele()) {
                        values[minor] += q;
                        countMinor++;
                    }
                }
            }
        }

        if (countMajor > 0)
            values[major] = values[major] / (double)countMajor;
        if (countMinor > 0)
            values[minor] = values[minor] / (double)countMinor;

        return values;
    }
}
