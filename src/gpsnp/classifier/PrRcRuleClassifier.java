package gpsnp.classifier;

import gpsnp.snpCalling.VariantCandidate;
import gpsnp.snpCalling.VariantConverter;
import snpsvm.bamreading.variant.Variant;

/**
 * Created by abrari on 27/07/15.
 */
public class PrRcRuleClassifier implements VariantClassifier {

    private boolean isTrueSNP(VariantCandidate v) {
        if (v.val("max.qual.minor") > 59.699 && v.val("total.depth") <= 82.149 && v.val("allele.balance") >= 0.067) {
            return true;
        } else if(v.val("max.qual.minor") <= 60.101 || v.val("total.depth") > 60.973 || (v.val("allele.balance") > 0.048 && v.val("allele.balance") < 0.930)) {
            return false;
        } else {
            return false;
        }
    }

    @Override
    public Variant classify(VariantCandidate v) {
        if (isTrueSNP(v)) {
            return VariantConverter.candidateToVariant(v);
        } else {
            return null;
        }
    }

}
