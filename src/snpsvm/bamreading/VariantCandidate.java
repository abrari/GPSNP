package snpsvm.bamreading;

import gpsnp.app.FeatureList;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.List;

/**
 * A class to store variant candidate (variant site) found during reading alignment column
 *
 * Created by abrari on 4/29/14.
 */
public class VariantCandidate {

    private int position;
    private char refBase;
    private FastaWindow refReader;
    private AlignmentColumn alnCol;
    private double[] featureValues;
    private int flankLeft;
    private int flankRight;

    private DecimalFormat formatter = new DecimalFormat("0.0####");

    public VariantCandidate(int position, char refBase, FastaWindow refReader, AlignmentColumn alnCol) {
        this.position = position;
        this.refBase = refBase;
        this.refReader = refReader;
        this.alnCol = alnCol;
        flankLeft = flankRight = 0;
    }

    public void computeFeatures() {
        List<FeatureComputer> computers = FeatureList.getFeatures();
        featureValues = new double[computers.size()];

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

    public void printFeatures(PrintStream out) {
        for (double value : featureValues) {
            out.print("\t" + formatter.format(value));
        }
    }

    public int getPosition() {
        return position;
    }


}
