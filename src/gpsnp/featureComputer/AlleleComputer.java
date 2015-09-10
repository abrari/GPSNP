package gpsnp.featureComputer;

import snpsvm.bamreading.FeatureComputer;

/**
 * Abstract class for allele-based computer (major and minor allele)
 *
 * Created by abrari on 4/28/14.
 */
public abstract class AlleleComputer implements FeatureComputer {

    protected boolean isPhred33 = false;

    final double[] values = new double[2];

    static final int major = 0;
    static final int minor = 1;

    @Override
    public int getColumnCount() {
        return values.length;
    }

    private char majorAllele;
    private char minorAllele;

    protected char getMajorAllele() {
        return majorAllele;
    }

    protected char getMinorAllele() {
        return minorAllele;
    }

    /**
     * Calculate major and minor alleles
     *
     * @param bases
     */
    protected void calculateAlleles(byte[] bases, int depth, char refBase) {

        // new concept: major allele = ref base, minor allele = alt base

        int[] baseCount = new int[4];
        baseCount[0] = 0; // A
        baseCount[1] = 0; // C
        baseCount[2] = 0; // G
        baseCount[3] = 0; // T

        for (int i=0; i < depth; i++) {
            byte base = bases[i];
            if ((char)base == 'A') baseCount[0]++;
            if ((char)base == 'C') baseCount[1]++;
            if ((char)base == 'G') baseCount[2]++;
            if ((char)base == 'T') baseCount[3]++;
        }

        majorAllele = refBase;
        minorAllele = computeAlt(refBase, baseCount);

    }

    private static char computeAlt(char ref, int[] counts) {
        int A = counts[0];
        int C = counts[1];
        int G = counts[2];
        int T = counts[3];

        if (ref == 'A')
            A = 0;
        if (ref == 'C')
            C = 0;
        if (ref == 'G')
            G = 0;
        if (ref == 'T')
            T = 0;

        //Find max of all...
        if (A >= C && A>=G && A>=T) {
            return 'A';
        }
        if (C >= A && C>=G && C>=T) {
            return 'C';
        }
        if (G >= A && G>=C && G>=T) {
            return 'G';
        }
        if (T >= A && T>=C && T>=G) {
            return 'T';
        }

        return 'N';
    }

}

