package gpsnp.classifier;

import gpsnp.snpCalling.VariantCandidate;
import snpsvm.bamreading.variant.Variant;

/**
 * Classification rule that maximizes sensitivity and specificity.
 * Rule was generated using Genetic Programming.
 *
 * Created by abrari on 26/07/15.
 */

public class SeSpRuleClassifier implements VariantClassifier {

    private boolean isTrueSNP(VariantCandidate v) {
        if (v.val("max.qual.minor") >= 58.341 && v.val("max.qual.major") >= 46.184 && v.val("total.depth") < 47.256) {
            return true;
        } else if(v.val("max.qual.minor") <= 58.933 || v.val("max.qual.major") < 48.921 || v.val("total.depth") >= 63.330) {
            return false;
        } else {
            return false;
        }
    }

    @Override
    public Variant classify(VariantCandidate v) {

        if (isTrueSNP(v)) {
            return new Variant(v.getContig(), v.getPosition(), v.getRefBase(), 'A', 0.0, 10, 2, 0.0, 0.0, 0.0);

        } else {
            return null;
        }
    }


}






















