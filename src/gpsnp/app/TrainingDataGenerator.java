package gpsnp.app;

import snpsvm.app.ArgParser;
import snpsvm.bamreading.FeatureComputer;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Generate training data (CSV-like)
 * Here we compute flank left and right, and determine training class
 *
 * @author abrari
 */
public class TrainingDataGenerator {

    private static String columnSeparator = ",";

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
        }
    }

    private static boolean candidateIsTrue(String contig, Integer pos) {
        if(trueSnps.containsKey(contig)) {
            return trueSnps.get(contig).contains(pos);
        }
        return false;
    }

    private static void printCsvHeader() {
        for(FeatureComputer feature : FeatureList.getFeatures()) {
            for(int i=0; i<feature.getColumnCount(); i++) {
                System.out.print(feature.getName(i) + columnSeparator);
            }
        }
        System.out.println("flank.left" + columnSeparator + "flank.right" + columnSeparator + "class");
    }

    private static void printArffHeader() {
        System.out.println("@RELATION gpsnp");

        for(FeatureComputer feature : FeatureList.getFeatures()) {
            for(int i=0; i<feature.getColumnCount(); i++) {
                if(feature.getName(i).equals("ts.tv"))
                    System.out.println("@ATTRIBUTE\t" + feature.getName(i) + "\t{ts,tv}");
                else
                    System.out.println("@ATTRIBUTE\t" + feature.getName(i) + "\tNUMERIC");
            }
        }

        System.out.println("@ATTRIBUTE\tflank.left\tNUMERIC");
        System.out.println("@ATTRIBUTE\tflank.right\tNUMERIC");
        System.out.println("@ATTRIBUTE\tclass\t{true,false}");

        System.out.println("@DATA");
    }

    public static void main(String[] args) throws Exception {

        ArgParser inputParser = new ArgParser(args);

        String trueSnpPath = inputParser.getStringArg("-T");
        readTrueSnps(trueSnpPath);

        String candidateSnpPath = inputParser.getStringArg("-C");   // Output from gpsnp.app.FeatureComputation
        BufferedReader reader = new BufferedReader(new FileReader(candidateSnpPath));
        String line;

        String headerStyle = inputParser.getStringArg("-h");

        // Reading all positions
        // Key: contig name, value: list of positions
        HashMap<String, ArrayList<Integer>> contigPos = new HashMap<String, ArrayList<Integer>>();
        HashMap<String, Integer> contigPosIdx = new HashMap<String, Integer>();

        while((line = reader.readLine()) != null) {
            String[] d = line.split("\t");
            if (d[0].equalsIgnoreCase("Chr")) continue;

            ArrayList<Integer> positions = contigPos.get(d[0]);
            if(positions == null) {
                positions = new ArrayList<Integer>();
                positions.add(Integer.valueOf(d[1]));
                contigPos.put(d[0], positions);
                contigPosIdx.put(d[0], new Integer(0));
            } else {
                positions.add(Integer.valueOf(d[1]));
            }
        }

        // Read again for feature values
        reader.close();
        reader = new BufferedReader(new FileReader(candidateSnpPath));

        // Print header
        if (headerStyle.equalsIgnoreCase("csv")) {
            printCsvHeader();
        } else if (headerStyle.equalsIgnoreCase("arff")) {
            printArffHeader();
        }

        while((line = reader.readLine()) != null) {
            String[] d = line.split("\t");
            String chr = d[0];
            if (chr.equalsIgnoreCase("Chr")) continue;

            // Flank size computation
            int pos = Integer.parseInt(d[1]);
            int curPosIdx = contigPosIdx.get(chr).intValue();
            int prevPos = contigPos.get(chr).get(Math.max(curPosIdx - 1, 0));
            int nextPos = contigPos.get(chr).get(Math.min(curPosIdx + 1, contigPos.get(chr).size() - 1));

            int flankLeft = pos - prevPos;
            int flankRight = nextPos - pos;

            // Echo-ing the features
            for (int i = 4; i < d.length - 2; i++) {
                if(i == 4) {
                    if(d[i].equals("0.0"))
                        System.out.print("ts");
                    else
                        System.out.print("tv");
                } else {
                    System.out.print(d[i]);
                }
                System.out.print(columnSeparator);
            }
            System.out.print(flankLeft + columnSeparator + flankRight + columnSeparator);

            if(candidateIsTrue(chr, pos)) {
                System.out.print("true");
            } else {
                System.out.print("false");
            }

            System.out.println();

            contigPosIdx.put(chr, curPosIdx+1);
        }

    }

}
