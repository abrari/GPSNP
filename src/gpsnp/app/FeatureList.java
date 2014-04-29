package gpsnp.app;

import gpsnp.featureComputer.*;
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
        computers.add(new MaxQualAlleleComputer());
        computers.add(new MeanQualAlleleComputer());
        computers.add(new FreqAlleleComputer());
        computers.add(new RelativeDistanceComputer());

        return computers;
    }

}
