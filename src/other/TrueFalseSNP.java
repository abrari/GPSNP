package other;

import java.io.BufferedReader;
import java.io.FileReader;
import java.text.DecimalFormat;
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

        // Key: contig name, value = [0: trueCount, 1: falseCount]
        HashMap<String, int[]> contigTFCount = new HashMap<String, int[]>();

        while((line = reader.readLine()) != null) {

            String[] d = line.split("\t");
            if(d[0].equalsIgnoreCase("Chr")) continue;

            System.out.print(d[0] + ":" + d[1] + "\t");

            if(candidateIsTrue(d[0], Integer.valueOf(d[1]))) {
                System.out.print("✓");

                int[] tf = contigTFCount.get(d[0]);
                if(tf == null) {
                    tf = new int[2];
                    tf[0]++;
                    contigTFCount.put(d[0], tf);
                } else {
                    tf[0]++;
                }
            } else {
                int[] tf = contigTFCount.get(d[0]);
                if(tf == null) {
                    tf = new int[2];
                    tf[1]++;
                    contigTFCount.put(d[0], tf);
                } else {
                    tf[1]++;
                }
            }

            System.out.println();

        }

        DecimalFormat df = new DecimalFormat("#.00");

        for (Map.Entry<String, int[]> entry : contigTFCount.entrySet()) {
            int trueCount = entry.getValue()[0];
            int falseCount = entry.getValue()[1];
            System.err.println("Statistics for " + entry.getKey() + ":");
            System.err.println("- True variants:  " + trueCount + " (" +  df.format((double)trueCount/(trueCount+falseCount)*100) + "%)");
            System.err.println("- False variants: " + falseCount + " (" + df.format((double)falseCount/(trueCount+falseCount)*100) + "%)");
        }

    }

}
