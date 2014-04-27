package gpsnp.app;

import snpsvm.bamreading.FeatureComputer;
import snpsvm.counters.DepthComputer;

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
        computers.add(new DepthComputer());

        return computers;
    }

}
