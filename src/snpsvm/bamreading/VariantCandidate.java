package snpsvm.bamreading;

import gpsnp.featureComputer.FeatureList;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.List;

/**
 * A class to store variant candidate (variant site) found during reading alignment column
 *
 * Created by abrari on 4/29/14.
 */
public class VariantCandidate {

    private String contig;
    private int position;
    private char refBase;
    private FastaWindow refReader;
    private AlignmentColumn alnCol;
    private double[] featureValues;

    private DecimalFormat formatter = new DecimalFormat("0.0####");

    public VariantCandidate(String contig, int position, char refBase, FastaWindow refReader, AlignmentColumn alnCol) {
        this.contig = contig;
        this.position = position;
        this.refBase = refBase;
        this.refReader = refReader;
        this.alnCol = alnCol;
    }

    public void computeFeatures(List<FeatureComputer> computers) {
        int featureCount = 0;
        for(FeatureComputer feature : computers) {
            for(int i=0; i<feature.getColumnCount(); i++) {
                featureCount++;
            }
        }

        featureValues = new double[featureCount];

        int c = 0;
        for(FeatureComputer computer: computers) {
            final double[] values = computer.computeValue(refBase, refReader, alnCol);
            for(int i=0; i<values.length; i++) {
                if (Double.isInfinite(values[i]) || Double.isNaN(values[i])) {
                    throw new IllegalArgumentException("Non-regular value for counter: " + computer.getName(i) + " found value=" + values[i]);
                }
                featureValues[c] = values[i];
                c++;
            }
        }
    }

    public void printMetaData(PrintStream out) {
        out.print(contig + "\t" + alnCol.getCurrentPosition() + "\t" + refBase);
    }

    public void printFeatureValues(PrintStream out) {
        for (double value : featureValues) {
            out.print("\t" + formatter.format(value));
        }
    }

    public int getPosition() {
        return position;
    }

    public char getRefBase() {
        return refBase;
    }

}
