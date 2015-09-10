package gpsnp.featureComputer;

import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FastaWindow;
import snpsvm.bamreading.MappedRead;

import java.util.Iterator;

/**
 * Maximum quality of major and minor allele
 *
 * Created by abrari on 4/28/14.
 */
public class MaxQualAlleleComputer extends AlleleComputer {

    public MaxQualAlleleComputer() {

    }

    public MaxQualAlleleComputer(boolean isPhred33) {
        this.isPhred33 = isPhred33;
    }

    @Override
    public String getName(int which) {
        if (which == major)
            return "max.qual.major";
        else
            return "max.qual.minor";
    }

    @Override
    public String getColumnDesc(int which) {
        if (which == major)
            return "Maximum quality of major allele";
        else
            return "Maximum quality of minor allele";
    }

    @Override
    public double[] computeValue(char refBase, FastaWindow window, AlignmentColumn col) {

        values[major] = 0.0;
        values[minor] = 0.0;

        int maxQualMajor = 0;
        int maxQualMinor = 0;

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
                        if (q > maxQualMajor) maxQualMajor = q;
                    } else if ((char) b == getMinorAllele()) {
                        if (q > maxQualMinor) maxQualMinor = q;
                    }
                }
            }
        }

        values[major] = maxQualMajor;
        values[minor] = maxQualMinor;

        return values;
    }
}
