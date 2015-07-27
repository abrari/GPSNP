package gpsnp.featureComputer;

import gpsnp.featureComputer.*;
import snpsvm.bamreading.FeatureComputer;
import snpsvm.counters.*;

import java.util.ArrayList;
import java.util.List;

/**
 * List of features to be computed
 *
 * Created by abrari on 4/27/14.
 */
public class FeatureList {

    public static List<FeatureComputer> getFeatures() {
        List<FeatureComputer> computers = new ArrayList<FeatureComputer>();
        computers.add(new TsTvComputer());
        computers.add(new MaxQualAlleleComputer());
        computers.add(new MeanQualAlleleComputer());
        computers.add(new RelativeDistanceComputer());
        computers.add(new DepthComputer());
        computers.add(new AlignmentQualityComputer());
        computers.add(new ErrorProbComputer());
        computers.add(new DinucRepeatCounter());
        computers.add(new StrandBiasComputer());
        computers.add(new AreaMismatchComputer());
        computers.add(new HomopolymerRunComputer());
        computers.add(new NucDiversityComputer());
        computers.add(new MismatchComputer());
        computers.add(new AlleleBalanceComputer());
        computers.add(new NearbyQualComputer());

        return computers;
    }

    // Exclude features that don't exists in the rule set
    public static List<FeatureComputer> getFeatures(boolean isPhred33) {
        List<FeatureComputer> computers = new ArrayList<FeatureComputer>();
        // computers.add(new TsTvComputer());
        computers.add(new MaxQualAlleleComputer(isPhred33));
        computers.add(new MeanQualAlleleComputer(isPhred33));
        // computers.add(new RelativeDistanceComputer());
        computers.add(new DepthComputer());
        // computers.add(new AlignmentQualityComputer());
        // computers.add(new ErrorProbComputer());
        // computers.add(new DinucRepeatCounter());
        // computers.add(new StrandBiasComputer());
        // computers.add(new AreaMismatchComputer());
        // computers.add(new HomopolymerRunComputer());
        // computers.add(new NucDiversityComputer());
        // computers.add(new MismatchComputer());
        computers.add(new AlleleBalanceComputer());
        // computers.add(new NearbyQualComputer(isPhred33));

        return computers;
    }

    public static int getFeatureCount() {
        int c = 0;
        for(FeatureComputer feature : getFeatures()) {
            for(int i=0; i<feature.getColumnCount(); i++) {
                c++;
            }
        }
        return c;
    }

    public static void printNames() {
        printNames("\t");
    }

    public static void printNames(String columnSeparator) {
        for(FeatureComputer feature : getFeatures()) {
            for(int i=0; i<feature.getColumnCount(); i++) {
                System.out.print(columnSeparator + feature.getName(i));
            }
        }
    }

    public static void main(String[] args) {
        for(FeatureComputer feature : getFeatures()) {
            for(int i=0; i<feature.getColumnCount(); i++) {
                System.out.println(feature.getName(i) + "\t" + feature.getClass().getSimpleName());
            }
        }
    }

}
