package gpsnp.snpCalling;

import gpsnp.featureComputer.FeatureList;
import snpsvm.bamreading.*;
import snpsvm.bamreading.intervalProcessing.IntervalCaller;
import snpsvm.bamreading.intervalProcessing.IntervalList;
import snpsvm.bamreading.intervalProcessing.IntervalList.Interval;
import snpsvm.bamreading.variant.Variant;

import java.io.*;
import java.util.List;

/**
 * Created by abrari on 26/07/15.
 */
public class GPSNPCaller implements IntervalCaller<List<Variant>> {

    protected static int instanceCount = 0;
    protected int myNumber = instanceCount;

    protected final File referenceFile;
    protected final IntervalList intervals;
    protected List<FeatureComputer> counters;
    protected List<Variant> variants = null;
    protected BAMWindowStore bamWindows;
    protected CallingOptions options = null;

    private long basesComputed = 0;

    public GPSNPCaller(File referenceFile,
                     IntervalList intervals,
                     BAMWindowStore bamWindows,
                     CallingOptions options) {
        this.referenceFile = referenceFile;
        this.intervals = intervals;
        this.counters = FeatureList.getFeatures();
        this.bamWindows = bamWindows;
        this.options = options;
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
            ReferenceBAMEmitter emitter = new ReferenceBAMEmitter(referenceFile, counters, window, options);

            //Store intermediate results in temporary files
            String tmpDataPrefix = "." + generateRandomString(12);

            File data = new File(tmpDataPrefix + ".tmp");

            //Read BAM file, write results to temporary file
            PrintStream dataStream = new PrintStream(new FileOutputStream(data));

            for(String contig : intervals.getContigs()) {
                for(Interval interval : intervals.getIntervalsInContig(contig)) {
                    emitter.emitWindow(contig, interval.getFirstPos(), interval.getLastPos(), dataStream);
                    basesComputed += interval.getSize();
                }
            }

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
