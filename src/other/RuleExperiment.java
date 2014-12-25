package other;

import com.opencsv.*;

import java.io.FileNotFoundException;
import java.io.FileReader;

public class RuleExperiment {

    // Feature indexes
    private static int idx_max_q_minor     = 2;
    private static int idx_mean_q_minor    = 4;
    private static int idx_total_depth     = 6;
    private static int idx_mapping_qual    = 7;
    private static int idx_error_prob      = 8;
    private static int idx_allele_balance  = 15;
    private static int idx_nuc_diversity   = 13;

    private static double[] mapToDouble(String[] strings) {
        double[] result = new double[strings.length-1]; // skip "class" attribute

        for(int i = 1; i < strings.length-1; i++) {
            result[i] = Double.parseDouble(strings[i]);
        }

        return result;
    }

    public static void main(String[] args) throws Exception {

        String arffFile = "/home/abrari/Tesis/GP-Experiment/data/Gm01.arff";

        CSVReader reader = new CSVReader(new FileReader(arffFile), ',', '\'', 21);
        String[] line;
        double[] d;

        // 0 = true, 1 = false
        int predicted;
        int actual;

        int i = 0;

        int[][] confMatrix = new int[2][2];

        while ((line = reader.readNext()) != null) {
            d = mapToDouble(line);
            actual = line[line.length-1].equals("true") ? 0 : 1;

            // RULE: FP2
            if (d[idx_allele_balance] >= 0.047743 && d[idx_total_depth] <= 123.902467 && d[idx_max_q_minor] > 61.171927) {
                predicted = 0;
            } else if (d[idx_nuc_diversity] < 0.508747 && d[idx_max_q_minor] < 60.277954) {
                predicted = 1;
            } else {
                predicted = 1;
            }

            confMatrix[actual][predicted]++;

            // if (i++ > 500) break;
        }

        System.out.println(confMatrix[0][0] + "\t" + confMatrix[0][1]);
        System.out.println(confMatrix[1][0] + "\t" + confMatrix[1][1]);
    }
}
