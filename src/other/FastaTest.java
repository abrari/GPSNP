package other;

import snpsvm.bamreading.FastaIndex;
import snpsvm.bamreading.FastaReader2;
import snpsvm.bamreading.FastaWindow;

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

        File fastaFile = new File("/home/abrari/Tesis/sam/lambda.fa");

        try {
            FastaReader2 fastaReader = new FastaReader2(fastaFile);

            Collection<String> contigs = fastaReader.getIndex().getContigs();

            System.out.println("Available contigs:");
            for (String contig : contigs) {
                System.out.println("\t" + contig + "\tlength: " + fastaReader.getContigLength(contig));
            }

            fastaReader.advanceToContig("gi");

            for (int i = 0; i < 10; i++) {
                fastaReader.advanceToPosition(i);
                System.out.print(fastaReader.nextBase());
            }

            System.out.println();

        } catch (Exception e) {
            e.printStackTrace();
        }


//        try {
//            FastaWindow window = new FastaWindow(fastaFile);
//
//            for (String contig : window.getContigs()) {
//                //System.out.println(contig);
//
//                // this will only read as long as FastaWindow::windowSize
//                window.resetTo(contig, 1);
//                System.out.println(window.allToString());
//            }
//
//        } catch (Exception e) {
//            e.printStackTrace();
//        }


    }

}
