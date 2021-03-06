package snpsvm.app;

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

import javax.swing.Timer;

import libsvm.LIBSVMModel;
import snpsvm.bamreading.BAMWindowStore;
import snpsvm.bamreading.CallingOptions;
import snpsvm.bamreading.FastaIndex.IndexNotFoundException;
import snpsvm.bamreading.FastaReader2;
import snpsvm.bamreading.FeatureComputer;
import snpsvm.bamreading.intervalProcessing.IntervalList;
import snpsvm.bamreading.intervalProcessing.IntervalList.Interval;
import snpsvm.bamreading.snpCalling.IntervalSNPCaller;
import snpsvm.bamreading.variant.VCFVariantEmitter;
import snpsvm.bamreading.variant.Variant;
import snpsvm.counters.CounterSource;

/**
 * Module that calls SNPs from an existing model. 
 * @author brendan
 *
 */
public class Predictor extends AbstractModule {

	private boolean emitProgress = true;
	
	@Override
	public boolean matchesModuleName(String name) {
		return name.equalsIgnoreCase("predict") || name.equals("emit");
	}
	
	public void emitColumnNames() {
		int index = 1;
		for(FeatureComputer counter : CounterSource.getCounters()) {
			for(int i=0; i<counter.getColumnCount(); i++) {
				System.out.println(index + "\t" + counter.getColumnDesc(i));
				index++;
			}
		}
	}

	@Override
	public void performOperation(String name, ArgParser args) {
		if (name.equals("emit")) {
			emitColumnNames();
			return;
		}
		
		if (CommandLineApp.configModule.getProperty("libsvm") == null) {
			System.out.println("\n  To begin, you must install libsvm. It's freely available from : http://www.csie.ntu.edu.tw/~cjlin/libsvm/");
			System.out.println("  Once you have downloaded and installed libsvm, tell SNPSVM where to find it, like this: ");
			System.out.println("  java snpsvm.jar config -add libsvm=/path/to/libsvm/ ");
			return;
		}
		
		String referencePath;
		String inputBAMPath;
		String modelPath;
		String vcfPath;
		
		try {
			referencePath = getRequiredStringArg(args, "-R", "Missing required argument for reference file, use -R");
			inputBAMPath = getRequiredStringArg(args, "-B", "Missing required argument for input BAM file, use -B");
			modelPath = getRequiredStringArg(args, "-M", "Missing required argument for model file, use -M");
			vcfPath = getRequiredStringArg(args, "-V", "Missing required argument for destination vcf file, use -V");
		} catch (MissingArgumentException e1) {
			System.err.println(e1.getMessage());
			return;
		}
		
		//Mostly for debugging, allows user-specified exclusion of counters
		super.processExcludedIntervals(args);
		
		boolean writeData = ! args.hasOption("-X");
		if (!writeData) {
			System.err.println("Skipping reading of BAM file... re-calling variants from existing output");
		}
		
		
		IntervalList intervals = getIntervals(args);
		
		
		
		File inputBAM = new File(inputBAMPath);
		File reference = new File(referencePath);
		File model = new File(modelPath);
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
		if (!model.exists()) {
			System.err.println("Model file " + model.getAbsolutePath() + " not found");
			return;
		}

		
		//If no interval list supplied create one from the reference
		if (intervals == null) {
			try {
				intervals = (new FastaReader2(reference)).toIntervals();
			} catch (IOException e) {
				System.err.println("There was an error reading the reference file, cannot proceed.");
			} catch (IndexNotFoundException e) {
				System.err.println("No index found for the reference file, cannot proceed.");
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
		
		emitProgress = ! args.hasOption("-quiet");
		
		ops.setRemoveTempFiles( ! args.hasOption("-preserve") );
		
		try {
			callSNPs(inputBAM, reference, model, vcf, intervals, ops);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.err.println("There was an error reading some of the files, cannot proceed.");
		} catch (IndexNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.err.println("Could not find the index file for the reference, please create one using samtools or picard");
		}
	}
	
	/**
	 * Create a new interval list with no intervals extending beyond the range given
	 * in the reference file
	 * @param reference
	 * @param intervals
	 * @return
	 * @throws IOException 
	 * @throws IndexNotFoundException 
	 */
	protected static IntervalList validateIntervals(File reference, IntervalList intervals) throws IOException, IndexNotFoundException {
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
			for(Interval interval : intervals.getIntervalsInContig(contig)) {
				int newPos = (int) Math.min(maxLength, interval.getLastPos());
				newIntervals.addInterval(contig, interval.getFirstPos(), newPos);
			}
		}
		
		
		return newIntervals;
	}
	
	public void callSNPs(File inputBAM, 
			File ref,
			File model,
			File destination,
			IntervalList intervals,
			CallingOptions ops) throws IOException, IndexNotFoundException {
		
		if (intervals != null)
			intervals = validateIntervals(ref, intervals);
		
		int threads= CommandLineApp.configModule.getThreadCount();
		//Initialize BAMWindow store
		BAMWindowStore bamWindows = new BAMWindowStore(inputBAM, threads);
		
		Timer progressTimer = null;

		final long intervalExtent = intervals.getExtent();
		
		List<Variant> allVars; 

		ThreadPoolExecutor threadPool = (ThreadPoolExecutor) Executors.newFixedThreadPool(threads);

		//final SplitSNPAndCall caller = new SplitSNPAndCall(ref, bamWindows, model, threadPool, ops);
		final IntervalSNPCaller caller = new IntervalSNPCaller(threadPool, ops, ref, model, bamWindows);

		//Submit multiple jobs to thread pool, returns immediately
		caller.submitAll(intervals);

		if (emitProgress) {
			System.out.println("Calling SNPs over " + intervals.getExtent() + " bases with " + threads + " threads in " + caller.getCallerCount() + " chunks");

			progressTimer = new javax.swing.Timer(100, new ActionListener() {

				@Override
				public void actionPerformed(ActionEvent arg0) {
					emitProgressString(caller, intervalExtent);
				}
			});
			progressTimer.setDelay(419);
			progressTimer.start();
		}

		//Blocks until all variants are called
		allVars = caller.getResult();

		//Emit one more progress message
		if (emitProgress) {
			emitProgressString(caller, intervalExtent);
		}
		
		
		if (progressTimer != null)
			progressTimer.stop();
		
		Collections.sort(allVars);
		
		//Write the variants to a file
		PrintStream writer = new PrintStream(new FileOutputStream(destination));
		
		VCFVariantEmitter vcfWriter = new VCFVariantEmitter();
		try {
			vcfWriter.writeHeader(writer, new FastaReader2(ref), inputBAM.getName().replace(".bam", ""));
			vcfWriter.writeVariants(allVars, writer);
		} catch (IndexNotFoundException e) {
			e.printStackTrace();
		}
		
		writer.close();
	}

	
	@Override
	public void emitUsage() {
		System.out.println("Predictor (SNP caller) module");
		System.out.println(" -R reference file");
		System.out.println(" -B input BAM file");
		System.out.println(" -V output variant file");
		System.out.println(" -M model file produced by buildmodel");
		System.out.println(" ---- Optional arguments -----");
		System.out.println(" -q [1.0] minimum Phred-scaled quality to report variant");
		System.out.println(" -d [2] minimum total depth to examine for variant");
		System.out.println(" -v [2] minimum reads with variant allele required for variant calling");
		System.out.println(" -quiet [false] do not emit progress to std. out");
	}

}
