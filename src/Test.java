/**
 *
 * Testing samtools SDK
 *
 * Created by abrari on 4/2/14.
 */

import net.sf.samtools.*;

import java.io.File;

public class Test {

    public static void main(String args[])
    {
        File inputFile = new File("/home/abrari/Tesis/sam/output.bam");

        SAMFileReader input = new SAMFileReader(inputFile, new File(inputFile.getAbsolutePath() + ".bai"));
        AbstractBAMFileIndex bamIndex = (AbstractBAMFileIndex) input.getIndex();

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

    }

}
