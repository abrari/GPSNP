package gpsnp.snpCalling;

import snpsvm.bamreading.variant.Variant;
import snpsvm.util.BinomMath;

/**
 * Created by abrari on 27/07/15.
 */
public class VariantConverter {

    public static Variant candidateToVariant(VariantCandidate v) {

        char[] baseSymbol = {'A','C','G','T','N'};
        byte[] alnBases = v.getAlnBases();
        int[] baseCounts = countBases(alnBases);
        int altIndex = computeAlt(v.getRefBase(), baseCounts);
        char alt = baseSymbol[altIndex];
        int depth = baseCounts[0] + baseCounts[1] + baseCounts[2] + baseCounts[3];
        int varDepth = baseCounts[altIndex];

        int T = baseCounts[0] + baseCounts[1] + baseCounts[2] + baseCounts[3];
        int X = baseCounts[altIndex];

        if (T > 250) {
            X = (250 * X)/T;
            T = 250;
        }

        //Compute het prob
        //Each read has 50% chance of coming from source with a non-reference base
        double hetProb = BinomMath.binomPDF((int) Math.round(X), (int) Math.round(T), 0.5);

        //Compute homo non-reference prob
        double homNonRefProb = BinomMath.binomPDF((int) Math.round(X), (int) Math.round(T), 0.99);

        //Compute homo-reference prob
        double homRefProb = BinomMath.binomPDF((int) Math.round(X), (int) Math.round(T), 0.005);

        double qual = v.val("mean.qual.minor"); // Infer variant quality from mean qual of minor allele

        return new Variant(v.getContig(), v.getPosition(), v.getRefBase(), alt, qual, depth, varDepth, homRefProb, hetProb, homNonRefProb);

    }

    private static int[] countBases(byte[] alnBases) {
        int[] baseCount = new int[4];
        baseCount[0] = 0; // A
        baseCount[1] = 0; // C
        baseCount[2] = 0; // G
        baseCount[3] = 0; // T

        for (int i=0; i < alnBases.length; i++) {
            byte base = alnBases[i];
            if ((char)base == 'A') baseCount[0]++;
            if ((char)base == 'C') baseCount[1]++;
            if ((char)base == 'G') baseCount[2]++;
            if ((char)base == 'T') baseCount[3]++;
        }
        return baseCount;
    }

    private static int computeAlt(char ref, int[] counts) {
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
            return 0;
        }
        if (C >= A && C>=G && C>=T) {
            return 1;
        }
        if (G >= A && G>=C && G>=T) {
            return 2;
        }
        if (T >= A && T>=C && T>=G) {
            return 3;
        }

        return 4;
    }

}
