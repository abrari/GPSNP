package other;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

public class TrueFalseSNP {

    // Key: contig name, values: set of positions
    private static HashMap<String, HashSet<Integer>> trueSnps = new HashMap<String, HashSet<Integer>>();

    private static void readTrueSnps(String filepath) throws Exception {

        BufferedReader reader = new BufferedReader(new FileReader(filepath));
        String line;

        while((line = reader.readLine()) != null) {
            String[] lineSplit = line.split("\t");

            HashSet<Integer> positions = trueSnps.get(lineSplit[0]);

            if(positions == null) {
                positions = new HashSet<Integer>();
                positions.add(Integer.valueOf(lineSplit[1]));
                trueSnps.put(lineSplit[0], positions);
            } else {
                positions.add(Integer.valueOf(lineSplit[1]));
            }

            // trueSnps.put(lineSplit[0] + ":" + lineSplit[1], lineSplit[2] + "\t" + lineSplit[3]);
        }

    }

    private static void printTrueSnps() {
        for(Map.Entry<String, HashSet<Integer>> chrPos : trueSnps.entrySet()) {
            System.out.println(chrPos.getKey());
            for (Integer pos : chrPos.getValue()) {
                System.out.println(" - " + pos);
            }
        }
    }

    private static boolean candidateIsTrue(String contig, Integer pos) {

        if(trueSnps.containsKey(contig)) {
            if(trueSnps.get(contig).contains(pos)) {
                return true;
            } else {
                return false;
            }
        }

        return false;
    }

    public static void main(String[] args) throws Exception {

        String trueSnpPath = args[0];
        readTrueSnps(trueSnpPath);

        String candidateSnpPath = args[1];
        BufferedReader reader = new BufferedReader(new FileReader(candidateSnpPath));
        String line;

        int trueCount = 0;
        int falseCount = 0;

        while((line = reader.readLine()) != null) {

            String[] d = line.split("\t");
            System.out.print(d[0] + ":" + d[1] + "\t");

            if(candidateIsTrue(d[0], Integer.valueOf(d[1]))) {
                trueCount++;
                System.out.print("✓");
            } else {
                falseCount++;
            }

            System.out.println();

        }

        System.err.println("True variants:  " + trueCount + " (" + (double)trueCount/(trueCount+falseCount)*100 + "%)");
        System.err.println("False variants: " + falseCount + " (" + (double)falseCount/(trueCount+falseCount)*100 + "%)");

    }

}
