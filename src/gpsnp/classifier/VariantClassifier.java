package gpsnp.classifier;

import gpsnp.snpCalling.VariantCandidate;
import snpsvm.bamreading.variant.Variant;

/**
 * Classify a VariantCandidate using some GP-generated rules
 *
 * Created by abrari on 26/07/15.
 */
public interface VariantClassifier {

    public Variant classify(VariantCandidate variantCandidate);

}
