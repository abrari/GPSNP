package gpsnp.featureComputer;

import snpsvm.bamreading.FeatureComputer;

/**
 * Abstract class for allele-based computer (major and minor allele)
 *
 * Created by abrari on 4/28/14.
 */
public abstract class AlleleComputer implements FeatureComputer {

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
    protected void calculateAlleles(byte[] bases, int depth) {

        majorAllele = minorAllele = 'X';

        char[] baseSymbol = {'A','C','G','T'};
        int majorIndex, minorIndex, majorVal, minorVal;

        int[] count = new int[4];
        count[0] = 0; // A
        count[1] = 0; // C
        count[2] = 0; // G
        count[3] = 0; // T

        for (int i=0; i<depth; i++) {   // limit to depth, because AlignmentColumn.bases is overwritten in every run
            byte base = bases[i];
            if ((char)base == 'A') count[0]++;
            if ((char)base == 'C') count[1]++;
            if ((char)base == 'G') count[2]++;
            if ((char)base == 'T') count[3]++;
        }

        if (count[0] > count[1]) {
            minorVal = count[1];
            minorIndex = 1;
            majorVal = count[0];
            majorIndex = 0;
        } else {
            minorVal = count[0];
            minorIndex = 0;
            majorVal = count[1];
            majorIndex = 1;
        }

        for(int i = 2; i < count.length; i++){
            if(count[i] >= majorVal){
                minorVal=majorVal;
                minorIndex=majorIndex;
                majorVal=count[i];
                majorIndex=i;
            }
            else if(count[i] > minorVal){
                minorVal=count[i];
                minorIndex=i;
            }
        }

        majorAllele = baseSymbol[majorIndex];
        minorAllele = baseSymbol[minorIndex];

    }

}

