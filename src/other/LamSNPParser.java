package other;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
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

        String filepath = "/media/abrari/data2/soy/snps/Test.snp";

        BufferedReader reader = new BufferedReader(new FileReader(filepath));
        String line;
        String lineReplace;

        initNucleotides();

        while((line = reader.readLine()) != null) {
            lineReplace = line;
            for(Entry<String, String> nuc : multiNuc.entrySet()) {
                lineReplace = lineReplace.replace(nuc.getKey(), nuc.getValue());
            }
            System.out.println(lineReplace);
        }

        reader.close();

    }

}
