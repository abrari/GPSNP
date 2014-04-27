package gpsnp.app;

import java.io.File;
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

    public static void main(String[] args) throws Exception {

        ArgParser inputParser = new ArgParser(args);
        String referencePath = null;
        String inputBAMPath = null;
        IntervalList intervals;

        referencePath = inputParser.getStringArg("-R");
        inputBAMPath = inputParser.getStringArg("-B");

        // System.out.println("Using reference sequence :" + referencePath);
        // System.out.println("Using reads sequence     :" + inputBAMPath);

        File reference = new File(referencePath);
        File inputBAM = new File(inputBAMPath);

        intervals = (new FastaReader2(reference)).toIntervals();
        BamWindow window = new BamWindow(inputBAM);
        List<FeatureComputer> features = FeatureList.getFeatures();

        // Debug message

        System.out.print("Pos\tRef\tAlt");
        for(FeatureComputer feature : features) {
            for(int i=0; i<feature.getColumnCount(); i++) {
                System.out.print("\t" + feature.getName());
            }
        }
        System.out.println();

        ReferenceBAMEmitter emitter = new ReferenceBAMEmitter(reference, features, window, new CallingOptions());
        for(String contig : intervals.getContigs()) {
            for(IntervalList.Interval inter : intervals.getIntervalsInContig(contig)) {
                emitter.emitWindow(contig, inter.getFirstPos(), inter.getLastPos(), System.out);
            }
        }


    }

}
