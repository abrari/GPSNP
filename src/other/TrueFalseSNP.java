package other;

import java.io.BufferedReader;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import snpsvm.app.ArgParser;

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

        ArgParser inputParser = new ArgParser(args);

        String trueSnpPath = inputParser.getStringArg("-T");
        readTrueSnps(trueSnpPath);

        Integer printPositions = inputParser.getIntegerArg("-p");   // 0: just print statistics, 1: print all positions

        String candidateSnpPath = inputParser.getStringArg("-C");
        BufferedReader reader = new BufferedReader(new FileReader(candidateSnpPath));
        String line;

        // Key: contig name, value = [0: trueCount, 1: falseCount]
        HashMap<String, int[]> contigTFCount = new HashMap<String, int[]>();

        while((line = reader.readLine()) != null) {

            String[] d = line.split("\t");
            if(d[0].equalsIgnoreCase("Chr")) continue;

            if(printPositions == 1) System.out.print(d[0] + ":" + d[1] + "\t");

            if(candidateIsTrue(d[0], Integer.valueOf(d[1]))) {
                if(printPositions == 1) System.out.print("✓");

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

            if(printPositions == 1) System.out.println();

        }

        DecimalFormat df = new DecimalFormat("#.00");

        int totalTrueCount = 0, totalFalseCount = 0;

        for (Map.Entry<String, int[]> entry : contigTFCount.entrySet()) {
            int trueCount = entry.getValue()[0];
            int falseCount = entry.getValue()[1];
            totalTrueCount += trueCount;
            totalFalseCount += falseCount;
            System.err.println("Statistics for " + entry.getKey() + ":");
            System.err.println("- True variants:  " + trueCount + " (" +  df.format((double)trueCount/(trueCount+falseCount)*100) + "%)");
            System.err.println("- False variants: " + falseCount + " (" + df.format((double)falseCount/(trueCount+falseCount)*100) + "%)");
        }

        System.err.println("Total:");
        System.err.println("- True variants:  " + totalTrueCount + " (" +  df.format((double)totalTrueCount/(totalTrueCount+totalFalseCount)*100) + "%)");
        System.err.println("- False variants: " + totalFalseCount + " (" + df.format((double)totalFalseCount/(totalTrueCount+totalFalseCount)*100) + "%)");

    }

}
