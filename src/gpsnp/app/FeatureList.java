package gpsnp.app;

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
        computers.add(new ErrorProbComputer());
        computers.add(new DinucRepeatCounter());
        computers.add(new PosDevComputer());
        computers.add(new ReadPosCounter());

        return computers;
    }

}
