package gpsnp.app;

import gpsnp.featureComputer.BaseQualityComputer;
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
        computers.add(new BaseQualityComputer());
        computers.add(new NearbyQualComputer());

        return computers;
    }

    public static void printNames() {
        for(FeatureComputer feature : getFeatures()) {
            for(int i=0; i<feature.getColumnCount(); i++) {
                System.out.print("\t" + feature.getName(i));
            }
        }
        System.out.println("\tflank.left\tflank.right");
    }

}
