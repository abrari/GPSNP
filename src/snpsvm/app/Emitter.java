package snpsvm.app;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import snpsvm.bamreading.*;
import snpsvm.bamreading.FastaIndex.IndexNotFoundException;
import snpsvm.bamreading.intervalProcessing.IntervalList;
import snpsvm.bamreading.intervalProcessing.IntervalList.Interval;
import snpsvm.bamreading.FeatureComputer;
import snpsvm.counters.CounterSource;

public class Emitter extends AbstractModule {

	@Override
	public void emitUsage() {
		System.out.println(" Emitter: emit training data to terminal (for debugging only) ");
		System.out.println(" -R : Reference File ");
		System.out.println(" -B : BAM File ");
		System.out.println(" -L : intervals ");
	}

	@Override
	public boolean matchesModuleName(String name) {
		return name.startsWith("write");
	}

	@Override
	public void performOperation(String name, ArgParser args) {
		String referencePath = null;
		String inputBAMPath = null;
		try {
			referencePath = getRequiredStringArg(args, "-R", "Missing required argument for reference file, use -R");
			inputBAMPath = getRequiredStringArg(args, "-B", "Missing required argument for input BAM file, use -B");
		} catch (MissingArgumentException e1) {
			System.err.println(e1.getMessage());
			return;
		}
		
		IntervalList intervals = getIntervals(args);

		File referenceFile = new File(referencePath);
		File inputBAM = new File(inputBAMPath);
		List<FeatureComputer> counters = CounterSource.getCounters();
		BamWindow window = new BamWindow(inputBAM);
		CallingOptions ops = new CallingOptions();

        //If no interval list supplied create one from the reference
        if (intervals == null) {
            try {
                intervals = (new FastaReader2(referenceFile)).toIntervals();
            } catch (IOException e) {
                System.err.println("There was an error reading the reference file, cannot proceed.");
            } catch (IndexNotFoundException e) {
                System.err.println("No index found for the reference file, cannot proceed.");
            }
        }

        try {
			ReferenceBAMEmitter emitter = new ReferenceBAMEmitter(referenceFile, counters, window, ops);
			BufferedWriter posWriter = new BufferedWriter(new FileWriter("emitter.positions.txt"));
			emitter.setPositionsWriter( posWriter );
			for(String contig : intervals.getContigs()) {
                for(Interval inter : intervals.getIntervalsInContig(contig)) {
                    emitter.emitWindow(contig, inter.getFirstPos(), inter.getLastPos(), System.out);
				}
			}
			
			posWriter.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
            System.err.println("I/O Error");
            e.printStackTrace();
		} catch (IndexNotFoundException e) {
			// TODO Auto-generated catch block
            System.err.println("Index not found!");
			e.printStackTrace();
		}
		
	}

}
