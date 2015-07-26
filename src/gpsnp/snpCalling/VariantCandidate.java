package gpsnp.snpCalling;

import gpsnp.featureComputer.FeatureList;
import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.FastaWindow;
import snpsvm.bamreading.FeatureComputer;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
    private String alnBases;
    private Map<String, Double> featureValues; // key: feature name, value: feature values

    private DecimalFormat formatter = new DecimalFormat("0.0####");

    public VariantCandidate(String contig, int position, char refBase, FastaWindow refReader, AlignmentColumn alnCol) {
        this.contig = contig;
        this.position = position;
        this.refBase = refBase;
        this.refReader = refReader;
        this.alnCol = alnCol;
        this.alnBases = alnCol.getBasesAsString(); // .getBases doesn't work (all column identical), don't know why
        this.featureValues = new HashMap<String, Double>();
    }

    public void computeFeatures(List<FeatureComputer> computers) {
        for(FeatureComputer computer: computers) {
            final double[] values = computer.computeValue(refBase, refReader, alnCol);
            for(int i=0; i<values.length; i++) {
                if (Double.isInfinite(values[i]) || Double.isNaN(values[i])) {
                    throw new IllegalArgumentException("Non-regular value for counter: " + computer.getName(i) + " found value=" + values[i]);
                }
                featureValues.put(computer.getName(i) , values[i]);
            }
        }
    }

    public void printMetaData(PrintStream out) {
        out.print(contig + "\t" + alnCol.getCurrentPosition() + "\t" + refBase);
    }

    public void printFeatureValues(PrintStream out) {
        for (Double value : featureValues.values()) {
            out.print("\t" + formatter.format(value));
        }
    }

    public String getContig() {
        return contig;
    }

    public int getPosition() {
        return position;
    }

    public char getRefBase() {
        return refBase;
    }

    public Double val(String featureName) {
        Double val = this.featureValues.get(featureName);
        if (val != null)
            return val;
        else
            return 0.0;
    }

    public byte[] getAlnBases() {

        int depth = this.alnBases.length();
        byte[] bases = new byte[depth];

        for (int i = 0; i < depth; i++) {
            bases[i] = (byte)this.alnBases.charAt(i);
        }

        return bases;
    }

}
