package other;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

/**
 * Created by abrari on 19/05/14.
 */
public class LamSNPParser {

    private static HashMap<String, String> multiNuc;

    private static void initNucleotides() {

        multiNuc = new HashMap<String, String>();

        multiNuc.put("W", "A\tT");
        multiNuc.put("S", "C\tG");
        multiNuc.put("M", "A\tC");
        multiNuc.put("K", "G\tT");
        multiNuc.put("R", "A\tG");
        multiNuc.put("Y", "C\tT");

        multiNuc.put("B", "C\tG\tT");
        multiNuc.put("D", "A\tG\tT");
        multiNuc.put("H", "A\tC\tT");
        multiNuc.put("V", "A\tC\tG");

    }

    public static void main(String[] args) throws Exception {

        String filepath = args[0];

        BufferedReader reader = new BufferedReader(new FileReader(filepath));
        String line;

        initNucleotides();

        while((line = reader.readLine()) != null) {
            for(Entry<String, String> nuc : multiNuc.entrySet()) {
                line = line.replace(nuc.getKey(), nuc.getValue());  // replace multinucleotide
                line = line.replaceAll("-\t?", "");                 // replace deletions
            }
            String[] lineSplit = line.split("\t");

            String contig = lineSplit[0];
            String pos = lineSplit[1];
            String ref = lineSplit[2];
            String[] alts = Arrays.copyOfRange(lineSplit, 3, lineSplit.length);

            HashSet<String> tempSet = new HashSet<String>(Arrays.asList(alts));
            String[] uniqueAlts = tempSet.toArray(new String[tempSet.size()]);

            if(uniqueAlts.length == 0 || (uniqueAlts.length == 1 && ref.equals(uniqueAlts[0]))) {
                // skip where ref = alt
                continue;
            }

            System.out.print(contig + "\t" + pos + "\t" + ref + "\t");
            System.out.print(Arrays.toString(uniqueAlts));
            System.out.println();
        }

        reader.close();

    }

}
