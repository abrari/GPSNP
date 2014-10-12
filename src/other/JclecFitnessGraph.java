package other;

import snpsvm.app.ArgParser;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Created by abrari on 12/10/14.
 *
 * Generate graph data from JCLEC's fitness output
 *
 */
public class JclecFitnessGraph {

    // Key: generation, value: max fitness
    private static Map<Integer, Double> maxFitness = new LinkedHashMap<Integer, Double>();
    private static File reportDir;

    private static void readBojarczukReports() throws FileNotFoundException {
        List<Integer> generations = new ArrayList<Integer>();

        File[] reports = reportDir.listFiles(new FileFilter() {
            public boolean accept(File path) {
                return path.getName().startsWith("Iteration_0");
            }
        });

        for (File report : reports) {
            int gen = Integer.parseInt(report.getName().substring(report.getName().lastIndexOf("_") + 1).replaceAll(".rep", ""));
            generations.add(gen);
        }
        Collections.sort(generations);

        for (Integer generation : generations) {
            Scanner s = new Scanner(new File(reportDir + "/" + "Iteration_0_" + generation + ".rep"));
            double maxFitnessTrue = -1, maxFitnessFalse = -1;
            while (s.hasNext()) {
                String line = s.nextLine();
                if (line.indexOf("class = true") != -1 && maxFitnessTrue == -1) {
                    maxFitnessTrue = Double.parseDouble(line.substring(line.lastIndexOf("; Fitness: ") + 11));
                }
                if (line.indexOf("class = false") != -1 && maxFitnessFalse == -1) {
                    maxFitnessFalse = Double.parseDouble(line.substring(line.lastIndexOf("; Fitness: ") + 11));
                }
                if (maxFitnessTrue != -1 && maxFitnessFalse != -1) break;
            }

            double averageMaxFitness = (maxFitnessTrue + maxFitnessFalse) / 2.0;

            maxFitness.put(generation, averageMaxFitness);
        }
    }

    private static void readReports() throws FileNotFoundException {
        List<Integer> generations = new ArrayList<Integer>();

        File[] trueReports = reportDir.listFiles(new FileFilter() {
            public boolean accept(File path) {
                return path.getName().startsWith("Iteration_0");
            }
        });

        for (File trueReport : trueReports) {
            int gen = Integer.parseInt(trueReport.getName().substring(trueReport.getName().lastIndexOf("_") + 1).replaceAll(".rep", ""));
            generations.add(gen);
        }
        Collections.sort(generations);

        for (Integer generation : generations) {
            Scanner sTrue = new Scanner(new File(reportDir + "/" + "Iteration_0_" + generation + ".rep"));
            String firstLineTrue = sTrue.nextLine();
            double maxFitnessTrue = Double.parseDouble(firstLineTrue.substring(firstLineTrue.lastIndexOf("; Fitness: ") + 11));

            Scanner sFalse = new Scanner(new File(reportDir + "/" + "Iteration_1_" + generation + ".rep"));
            String firstLineFalse = sFalse.nextLine();
            double maxFitnessFalse = Double.parseDouble(firstLineFalse.substring(firstLineFalse.lastIndexOf("; Fitness: ") + 11));

            double averageMaxFitness = (maxFitnessTrue + maxFitnessFalse) / 2.0;

            maxFitness.put(generation, averageMaxFitness);
        }
    }

    public static void main(String[] args) throws FileNotFoundException {

        ArgParser inputParser = new ArgParser(args);
        boolean isBojarczuk = false;

        reportDir = new File(inputParser.getStringArg("-D"));
        File[] checkFiles = reportDir.listFiles(new FileFilter() {
            public boolean accept(File path) {
                return path.getName().startsWith("Iteration_1");
            }
        });

        if (checkFiles.length == 0) isBojarczuk = true;

        if (isBojarczuk) {
            readBojarczukReports();
        } else {
            readReports();
        }

        for (Map.Entry<Integer, Double> entry : maxFitness.entrySet()) {
            System.out.println(entry.getKey() + "\t" + entry.getValue());
        }

    }

}
