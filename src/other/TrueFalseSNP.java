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

    public static void main(String[] args) throws Exception {

        String trueSnpPath = args[0];
        readTrueSnps(trueSnpPath);

        for(Map.Entry<String, HashSet<Integer>> chrPos : trueSnps.entrySet()) {
            System.out.println(chrPos.getKey());
            for (Integer pos : chrPos.getValue()) {
                System.out.println(" - " + pos);
            }
        }

    }

}
