import snpsvm.bamreading.FastaIndex;
import snpsvm.bamreading.FastaReader2;

import java.io.File;
import java.util.Collection;

/**
 *
 * Testing reading FASTA files
 *
 * Created by abrari on 4/3/14.
 */
public class FastaTest {

    public static void main(String[] args) {

        File fastaFile = new File("/home/abrari/Tesis/sam/lambda_virus.fa");

        try {
            FastaReader2 fastaReader = new FastaReader2(fastaFile);

            Collection<String> contigs = fastaReader.getIndex().getContigs();

            System.out.println("Available contigs:");
            for (String contig : contigs) {
                System.out.println("\t" + contig + "\tlength: " + fastaReader.getContigLength(contig));
            }

            fastaReader.advanceToContig("gi|9626243|ref|NC_001416.1|");

            for (int i = 0; i < 5; i++) {
                fastaReader.advanceToPosition(i);
                System.out.print(fastaReader.nextBase());
            }

            System.out.println();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}
