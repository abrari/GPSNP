package snpsvm.bamreading;

/**
 * A few user-settable options for variant calling
 * @author brendan
 *
 */
public class CallingOptions {
	
	static final int DEFAULT_MIN_DEPTH = 2;
	static final int DEFAULT_MIN_VAR_DEPTH = 2;
	static final double DEFAULT_MIN_QUALITY = 1.0;
	
	int minTotalDepth = DEFAULT_MIN_DEPTH;
	int minVariantDepth = DEFAULT_MIN_VAR_DEPTH;
	double minQuality = DEFAULT_MIN_QUALITY;
	
	boolean removeTempFiles = true;
    boolean phred33Qual = false;    // is the input quality in phred+33 scoring?
	
	public CallingOptions() {
		
	}

    public boolean isPhred33Qual() {
        return phred33Qual;
    }

    public void setPhred33Qual(boolean phred33Qual) {
        this.phred33Qual = phred33Qual;
    }

    public boolean isRemoveTempFiles() {
		return removeTempFiles;
	}



	public void setRemoveTempFiles(boolean removeTempFiles) {
		this.removeTempFiles = removeTempFiles;
	}



	public int getMinTotalDepth() {
		return minTotalDepth;
	}

	public void setMinTotalDepth(int minTotalDepth) {
		this.minTotalDepth = minTotalDepth;
	}

	public int getMinVariantDepth() {
		return minVariantDepth;
	}

	public void setMinVariantDepth(int minVariantDepth) {
		this.minVariantDepth = minVariantDepth;
	}

	public double getMinQuality() {
		return minQuality;
	}

	public void setMinQuality(double minQuality) {
		this.minQuality = minQuality;
	}
	
	
}
