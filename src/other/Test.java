package other; /**
 *
 * Testing samtools SDK
 *
 * Created by abrari on 4/2/14.
 */

import net.sf.samtools.*;
import snpsvm.bamreading.AlignmentColumn;
import snpsvm.bamreading.BamWindow;
import snpsvm.bamreading.MappedRead;

import java.io.File;
import java.util.Iterator;

public class Test {

    public static void main(String args[])
    {
        File inputFile = new File("/home/abrari/Tesis/sam/output.bam");

        SAMFileReader input = new SAMFileReader(inputFile, new File(inputFile.getAbsolutePath() + ".bai"));
        AbstractBAMFileIndex bamIndex = (AbstractBAMFileIndex) input.getIndex();

        // getting contig names from BAM

        SAMFileHeader header = input.getFileHeader();
        SAMSequenceDictionary dictionary = header.getSequenceDictionary();

        for(SAMSequenceRecord seq : dictionary.getSequences()) {
            System.out.println(seq.getSequenceName());
        }


        /*

        int recordCount = 0;
        int mappedRecordCount = 0;

        for (SAMRecord samRecord : input)
        {
            recordCount++;
        }

        System.out.println("Number of records: " + recordCount);

        for (int i = 0; i < bamIndex.getNumberOfReferences(); i++)
        {
            mappedRecordCount += bamIndex.getMetaData(i).getAlignedRecordCount();
        }

        System.out.println("Number of mapped records: " + mappedRecordCount);

        BamWindow window = new BamWindow(inputFile);

        window.advanceToContig("gi|9626243|ref|NC_001416.1|");

        System.out.println(window.getCurrentPosition());
        System.out.println(window.getCurrentContig());

        AlignmentColumn column = new AlignmentColumn(window);

        for (int i=0; i<=250; i++) {
            column.advance();
            System.out.println((i+1) + "\t" + column.getBasesAsString());
        }

        */

    }

}
