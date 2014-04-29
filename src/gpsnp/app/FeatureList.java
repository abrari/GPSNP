package gpsnp.app;

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
        computers.add(new StrandBiasComputer());
        computers.add(new AreaMismatchComputer());
        computers.add(new HomopolymerRunComputer());
        computers.add(new NucDiversityComputer());
        computers.add(new MismatchComputer());
        computers.add(new AlleleBalanceComputer());

        return computers;
    }

}
