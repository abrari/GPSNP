package gpsnp.snpCalling;

import gpsnp.classifier.PrRcRuleClassifier;
import gpsnp.classifier.SeSpRuleClassifier;
import gpsnp.classifier.VariantClassifier;
import gpsnp.featureComputer.FeatureList;
import snpsvm.bamreading.*;
import snpsvm.bamreading.intervalProcessing.IntervalCaller;
import snpsvm.bamreading.intervalProcessing.IntervalList;
import snpsvm.bamreading.intervalProcessing.IntervalList.Interval;
import snpsvm.bamreading.variant.Variant;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by abrari on 26/07/15.
 */
public class GPSNPCaller implements IntervalCaller<List<Variant>> {

    private static int instanceCount = 0;
    private int myNumber = instanceCount;

    private final File referenceFile;
    private final IntervalList intervals;
    private List<FeatureComputer> counters;
    private List<Variant> variants = null;
    private BAMWindowStore bamWindows;
    private CallingOptions options = null;

    private VariantClassifier classifier;

    private long basesComputed = 0;

    public GPSNPCaller(File referenceFile,
                     IntervalList intervals,
                     BAMWindowStore bamWindows,
                     CallingOptions options) {
        this.referenceFile = referenceFile;
        this.intervals = intervals;
        this.counters = FeatureList.getFeatures(options.isPhred33Qual());
        this.bamWindows = bamWindows;
        this.options = options;
        this.classifier = loadClassifier(options.getClassifierClass());
        this.variants = new ArrayList<Variant>();
        instanceCount++;
    }

    @Override
    public List<Variant> getResult() {
        return variants;
    }

    @Override
    public boolean isResultReady() {
        return variants != null;
    }

    @Override
    public long getBasesCalled() {
        return basesComputed;
    }

    @Override
    public void run() {
        try {
            BamWindow window = bamWindows.getWindow();
            VariantCandidateEmitter emitter = new VariantCandidateEmitter(referenceFile, counters, window, options);
            List<VariantCandidate> variantCandidates = new ArrayList<VariantCandidate>();

            for(String contig : intervals.getContigs()) {
                for(Interval interval : intervals.getIntervalsInContig(contig)) {
                    List<VariantCandidate> candidates = emitter.emitWindow(contig, interval.getFirstPos(), interval.getLastPos());
                    if (candidates.size() > 0) {
                        variantCandidates.addAll(candidates);
                    }
                    basesComputed += interval.getSize();
                }
            }

            // Call true SNPs
            for (VariantCandidate v: variantCandidates) {
                // dataStream.println(v.getContig() + ":" + v.getPosition() + "\t" + v.getRefBase() + "\t" + Arrays.toString(v.getAlnBases()));
                Variant var = this.classifier.classify(v);
                if (var != null) {
                    this.variants.add(var);
                }
            }


        } catch (IOException iox) {
            iox.printStackTrace();
        } catch (FastaIndex.IndexNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private VariantClassifier loadClassifier(String className) {

        String prefix = "gpsnp.classifier.";

        try {
            Object classifier = Class.forName(prefix + className).newInstance();
            return (VariantClassifier)classifier;

        } catch (InstantiationException e) {
            e.printStackTrace();
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }

        return new PrRcRuleClassifier();

    }

}
