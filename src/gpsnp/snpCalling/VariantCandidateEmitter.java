package gpsnp.snpCalling;

/**
 * Created by abrari on 26/07/15.
 */
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import snpsvm.bamreading.*;
import snpsvm.bamreading.FastaIndex.IndexNotFoundException;
import snpsvm.bamreading.FastaReader2.EndOfContigException;

public class VariantCandidateEmitter {

    final FastaWindow refReader;
    protected AlignmentColumn alnCol;
    private Map<String, Integer> contigMap;
    List<FeatureComputer> counters;
    protected BufferedWriter positionWriter = null;
    protected DecimalFormat formatter = new DecimalFormat("0.0####");
    protected final int minDepth;
    protected final int minVarDepth;

    List<VariantCandidate> variants;

    public VariantCandidateEmitter(File reference, List<FeatureComputer> counters, BamWindow window, CallingOptions ops) throws IOException, IndexNotFoundException {
        refReader = new FastaWindow(reference);
        contigMap = new HashMap<String, Integer>();
        for(String contig : refReader.getContigs()) {
            contigMap.put(contig, refReader.getContigLength(contig).intValue());
        }
        alnCol = new AlignmentColumn(window);
        this.minDepth = ops.getMinTotalDepth();
        this.minVarDepth = ops.getMinVariantDepth();
        this.counters = counters;
        this.variants = new ArrayList<VariantCandidate>();
    }

    public VariantCandidate emitLine(String contig) {

        if (alnCol.getApproxDepth() >= minDepth) {
            final char refBase = refReader.getBaseAt(alnCol.getCurrentPosition());
            if (refBase == 'N') {
                return null;
            }

            boolean differringBases = alnCol.hasXDifferingBases(refBase, minVarDepth);

            if (! differringBases) {
                return null;
            }

            VariantCandidate var = new VariantCandidate(contig, alnCol.getCurrentPosition(), refBase, refReader, alnCol);
            var.computeFeatures(this.counters);

            return var;
        }

        return null;
    }

    public List<VariantCandidate> emitWindow(String contig, int start, int end) throws Exception {
        if (! refReader.containsContig(contig)) {
            //throw new IllegalArgumentException("Reference does not have contig : " + contig);
            System.err.println("Warning, reference does not contain contig: " + contig + ".  Skipping it.");
            return null;
        }

        if (! alnCol.containContig(contig)) {
            System.err.println("Warning, alignment does not contain contig: " + contig + ".  Skipping it.");
            return null;
        }

        try {
            refReader.resetTo(contig, Math.max(1, start-refReader.windowSize/2));
            alnCol.advanceTo(contig, start);

            int curPos = start;
            while(curPos < end && alnCol.hasMoreReadsInCurrentContig()) {
                VariantCandidate var = emitLine(contig);
                if (var != null)
                    this.variants.add(var);

                if (refReader.indexOfLeftEdge()<(alnCol.getCurrentPosition()-refReader.windowSize/2)) {
                    try {
                        refReader.shift();
                    }
                    catch(EndOfContigException ex) {
                        //don't worry about it
                    }
                }
                alnCol.advance(1);
                curPos++;

                //Sanity check
                if (alnCol.getCurrentPosition() != (curPos)) {
                    System.err.println("Yikes, bam reader position is not equal to current position");
                }
            }

            System.out.println("[" + contig + ":" + start + "-" + end + "] yields " + variants.size() + " variant candidates");

            return this.variants;

        } catch (EndOfContigException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        return null;
    }
}
