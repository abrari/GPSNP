package gpsnp.app;

import snpsvm.app.ArgParser;
import snpsvm.bamreading.BAMWindowStore;
import snpsvm.bamreading.CallingOptions;
import snpsvm.bamreading.FastaIndex;
import snpsvm.bamreading.FastaReader2;
import snpsvm.bamreading.intervalProcessing.IntervalList;
import snpsvm.bamreading.snpCalling.IntervalSNPCaller;
import snpsvm.bamreading.variant.VCFVariantEmitter;
import snpsvm.bamreading.variant.Variant;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;

/**
 * Created by abrari on 25/07/15.
 */
public class Main {

    // Arguments processing

    public static String getRequiredStringArg(ArgParser args, String arg, String errorMessage) throws Exception {
        if (args.hasOption(arg)) {
            return args.getStringArg(arg);
        }
        else {
            throw new Exception(errorMessage);
        }
    }

    public static String getOptionalStringArg(ArgParser args, String arg) {
        if (args.hasOption(arg)) {
            return args.getStringArg(arg);
        }
        else {
            return null;
        }
    }

    public static Double getOptionalDoubleArg(ArgParser args, String arg) {
        if (args.hasOption(arg)) {
            return args.getDoubleArg(arg);
        }
        else {
            return null;
        }
    }

    public static IntervalList getIntervals(ArgParser args) {
        String intervalsStr = getOptionalStringArg(args, "-L");
        if (intervalsStr == null) {
            return null;
        }
        IntervalList intervals = new IntervalList();
        File testFile = new File(intervalsStr);
        if (testFile.exists()) {
            try {
                intervals.buildFromBEDFile(testFile);
            } catch (IOException e) {
                System.err.println("Error building interval list from BED file :  " + e.getMessage());
            }
        }
        else {
            try {
                intervals.buildFromString(intervalsStr);
            }
            catch (Exception ex) {
                System.err.println("Error parsing intervals from " + intervalsStr + " : " + ex.getMessage());
                return null;
            }
        }

        return intervals;
    }

    /**
     * Create a new interval list with no intervals extending beyond the range given
     * in the reference file
     * @param reference
     * @param intervals
     * @return
     * @throws IOException
     * @throws snpsvm.bamreading.FastaIndex.IndexNotFoundException
     */
    protected static IntervalList validateIntervals(File reference, IntervalList intervals) throws IOException, FastaIndex.IndexNotFoundException {
        FastaReader2 refReader = new FastaReader2(reference);
        IntervalList newIntervals = new IntervalList();
        for(String contig : intervals.getContigs()) {
            if (! refReader.containsContig(contig)) {
                throw new IllegalArgumentException("No contig '" + contig + "' found in reference.");
            }
            Long maxLength = refReader.getContigLength(contig);
            if (maxLength == null) {
                throw new IllegalArgumentException("Could not read length of contig " + contig + " in reference");
            }
            for(IntervalList.Interval interval : intervals.getIntervalsInContig(contig)) {
                int newPos = (int) Math.min(maxLength, interval.getLastPos());
                newIntervals.addInterval(contig, interval.getFirstPos(), newPos);
            }
        }


        return newIntervals;
    }

    // Main

    public static void callSNPs(File inputBAM,
                                File ref,
                                File destination,
                                IntervalList intervals,
                                CallingOptions ops,
                                int threads) throws IOException, FastaIndex.IndexNotFoundException {

        if (intervals != null)
            intervals = validateIntervals(ref, intervals);

        //Initialize BAMWindow store
        BAMWindowStore bamWindows = new BAMWindowStore(inputBAM, threads);

        final long intervalExtent = intervals.getExtent();

        System.out.println("Range " + intervals.toString());

        List<Variant> allVars;

        ThreadPoolExecutor threadPool = (ThreadPoolExecutor) Executors.newFixedThreadPool(threads);

        final IntervalSNPCaller caller = new IntervalSNPCaller(threadPool, ops, ref, null, bamWindows);

        //Submit multiple jobs to thread pool, returns immediately
        caller.submitAll(intervals);

        //Blocks until all variants are called
        allVars = caller.getResult();

        Collections.sort(allVars);

        //Write the variants to a file
        PrintStream writer = new PrintStream(new FileOutputStream(destination));

        VCFVariantEmitter vcfWriter = new VCFVariantEmitter();
        try {
            vcfWriter.writeHeader(writer, new FastaReader2(ref), inputBAM.getName());
            vcfWriter.writeVariants(allVars, writer);
        } catch (FastaIndex.IndexNotFoundException e) {
            e.printStackTrace();
        }

        writer.close();
    }

    public static void main(String[] argv) {

        ArgParser args = new ArgParser(argv);

        String referencePath;
        String inputBAMPath;
        String vcfPath;
        Integer threads = 4;

        try {
            referencePath = getRequiredStringArg(args, "-R", "Missing required argument for reference file, use -R");
            inputBAMPath = getRequiredStringArg(args, "-B", "Missing required argument for input BAM file, use -B");
            vcfPath = getRequiredStringArg(args, "-V", "Missing required argument for destination vcf file, use -V");
        } catch (Exception e1) {
            System.err.println(e1.getMessage());
            return;
        }

        if (args.getIntegerArg("-t") != null) {
            threads = args.getIntegerArg("-t");
        }

        IntervalList intervals = getIntervals(args);

        File inputBAM = new File(inputBAMPath);
        File reference = new File(referencePath);
        File vcf = new File(vcfPath);

        //Some error checking...make sure files exist
        if (!inputBAM.exists()) {
            System.err.println("Input .BAM file " + inputBAM.getAbsolutePath() + " not found");
            return;
        }
        if (!reference.exists()) {
            System.err.println("Reference file " + reference.getAbsolutePath() + " not found");
            return;
        }

        //If no interval list supplied create one from the reference
        if (intervals == null) {
            try {
                intervals = (new FastaReader2(reference)).toIntervals();
            } catch (IOException e) {
                System.err.println("There was an error reading the reference file, cannot proceed.");
                System.exit(500);
            } catch (FastaIndex.IndexNotFoundException e) {
                System.err.println("No index found for the reference file, cannot proceed.");
                System.exit(404);
            }
        }

        //Generate a CallingOptions object with some params for variant calling
        CallingOptions ops = new CallingOptions();
        Double qCutoff = getOptionalDoubleArg(args, "-q");
        if (qCutoff != null)
            ops.setMinQuality(qCutoff);
        Double minDepthDub = getOptionalDoubleArg(args, "-d");
        if (minDepthDub != null) {
            ops.setMinTotalDepth((int) Math.floor(minDepthDub));
        }
        Double minVarDepthDub = getOptionalDoubleArg(args, "-v");
        if (minVarDepthDub != null) {
            ops.setMinVariantDepth((int)Math.floor(minVarDepthDub));
        }

        ops.setRemoveTempFiles( ! args.hasOption("-preserve") );

        try {
            callSNPs(inputBAM, reference, vcf, intervals, ops, threads);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
            System.err.println("There was an error reading some of the files, cannot proceed.");
        } catch (FastaIndex.IndexNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
            System.err.println("Could not find the index file for the reference, please create one using samtools or picard");
        }
    }

}
