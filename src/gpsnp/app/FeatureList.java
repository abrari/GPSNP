package gpsnp.app;

import gpsnp.featureComputer.AlignmentQualityComputer;
import gpsnp.featureComputer.BaseQualityComputer;
import snpsvm.bamreading.FeatureComputer;
import snpsvm.counters.DepthComputer;
import snpsvm.counters.TsTvComputer;

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
        computers.add(new TsTvComputer());
        computers.add(new AlignmentQualityComputer());
        computers.add(new BaseQualityComputer());

        return computers;
    }

}
