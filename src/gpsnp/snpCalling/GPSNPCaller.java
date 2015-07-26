package gpsnp.snpCalling;

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
        this.classifier = new SeSpRuleClassifier();
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

            //Temporary files for debugging
            String tmpDataPrefix = "." + generateRandomString(12);

            File data = new File(tmpDataPrefix + ".tmp");

            //Read BAM file, write results to temporary file
            PrintStream dataStream = new PrintStream(new FileOutputStream(data));

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
                    dataStream.println(var.toString());
                }
            }

            //CRITICAL: must return its bamWindow to the BAMWindowStore
            bamWindows.returnToStore(window);
            dataStream.close();

            //Remove temporary files
            if (options.isRemoveTempFiles()) {
                data.delete();
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

    private static String generateRandomString(int length) {
        StringBuilder strB = new StringBuilder();
        while(strB.length() < length) {
            char c = chars.charAt( (int)(chars.length()*Math.random()) );
            strB.append(c);
        }
        return strB.toString();
    }

    private static final String chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabscdefhgijklmnopqrstuvwxyz1234567890";

}
