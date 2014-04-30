package gpsnp.app;

import java.io.File;
import java.io.IOException;
import java.util.List;

import snpsvm.app.ArgParser;
import snpsvm.bamreading.*;
import snpsvm.bamreading.intervalProcessing.IntervalList;
import snpsvm.counters.CounterSource;

/**
 * Test computing features listed on FeatureList
 *
 * Created by abrari on 4/27/14.
 */
public class FeatureComputation {

    private static int minDepth = 2;
    private static int minVarDepth = 2;

    public static void main(String[] args) throws Exception {

        ArgParser inputParser = new ArgParser(args);
        String referencePath = null;
        String inputBAMPath = null;
        IntervalList intervals;

        referencePath = inputParser.getStringArg("-R");
        inputBAMPath = inputParser.getStringArg("-B");

        File reference = new File(referencePath);
        File inputBAM = new File(inputBAMPath);

        intervals = (new FastaReader2(reference)).toIntervals();

        FastaWindow refReader = new FastaWindow(reference);
        BamWindow window = new BamWindow(inputBAM);
        AlignmentColumn alnCol = new AlignmentColumn(window);

        List<FeatureComputer> features = FeatureList.getFeatures();

        // Debug message

        System.out.print("Pos\tRef\tAlt");
        FeatureList.printNames();

        for(String contig : intervals.getContigs()) {
            for(IntervalList.Interval inter : intervals.getIntervalsInContig(contig)) {
                int start = inter.getFirstPos();
                int end = inter.getLastPos();

                try {
                    refReader.resetTo(contig, Math.max(1, start - refReader.getWindowSize() / 2));
                    alnCol.advanceTo(contig, start);

                    int curPos = start;

                    VariantCandidate var = null;

                    int varCount = 0;
                    while(curPos < end && alnCol.hasMoreReadsInCurrentContig()) {
                        if (alnCol.getApproxDepth() >= minDepth) {
                            final char refBase = refReader.getBaseAt(alnCol.getCurrentPosition());
                            if (refBase != 'N' && alnCol.hasXDifferingBases(refBase, minVarDepth)) {
                                var = new VariantCandidate(alnCol.getCurrentPosition(), refBase, refReader, alnCol);

                                var.computeFeatures();
                                var.printMetaData(System.out);
                                var.printFeatureValues(System.out);

                                System.out.println();

                                varCount++;
                            }
                        }

                        if (refReader.indexOfLeftEdge()<(alnCol.getCurrentPosition()-refReader.getWindowSize()/2)) {
                            try {
                                refReader.shift();
                            }
                            catch(FastaReader2.EndOfContigException ex) {
                                //don't worry about it
                            }
                        }
                        alnCol.advance(1);
                        curPos++;
                    }
                } catch (FastaReader2.EndOfContigException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }


            }
        }


    }

}
