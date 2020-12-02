package uniandes.algorithms.tr;

import ngsep.genome.GenomicRegion;

public class TandemRepeat implements GenomicRegion{

	// first index of the tuple that has a match in last-d
	private int patternBeginning;

	// last index that indicates a match
	private int last;

	//length of the pattern
	private int unitLength;

	//sum of heads criteria for a tandem repeat
	private int sumOfHeads;

	//number of copies of the tandem repeat
	private double numCopies;

	//first index of the  tandem repeat
	private int first;

	//sequence name, provided by fasta file
	private String sequenceName;
	
	//quality score for diagnostic
	private short qualityScore;
	
	//pattern
	private String pattern;

	public TandemRepeat(int first, int last, int patternLength, int sumOfHeads) {
		this.first = first;
		this.last = last;
		this.unitLength = patternLength;
		this.sumOfHeads = sumOfHeads;
		this.numCopies = 0;
		this.patternBeginning = last - patternLength + 1;
		this.sequenceName = "";
		this.qualityScore = 0;
		this.pattern = "";
	}

	/**
	 * @return first index of the tuple that has a match in last-d
	 */
	public int getPatternBeginning() {
		return patternBeginning;
	}

	/**
	 * @return the last
	 */
	public int getLast() {
		return last;
	}

	/**
	 * @return the patternLength
	 */
	public int getUnitLength() {
		return unitLength;
	}

	/**
	 * @return the sumOfHeads
	 */
	public int getSumOfHeads() {
		return sumOfHeads;
	}

	/**
	 * @return the times
	 */
	public double getNumCopies() {
		return numCopies;
	}

	/**
	 * @param last the last to set
	 */
	public void setLast(int last) {
		this.last = last;
	}

	/**
	 * @param numCopies
	 *            the times to set
	 */
	public void setNumCopies(double numCopies) {
		this.numCopies = numCopies;
	}

	/**
	 * @return the beginning
	 */
	public int getFirst() {
		return first;
	}

	/**
	 * 
	 * @return distance used in TR candidate selector algorithm
	 */
	public int getDistance() {
		return last - patternBeginning + 1;
	}

	/**
	 * @param first
	 *            the beginning to set
	 */
	public void setFirst(int first) {
		if (first < 0) {
			this.first = 0;
		} else {
			this.first = first;
		}

	}

	/**
	 * 
	 * @return Total size of the TR
	 */
	public int getTotalSize() {
		return last - first + 1;
	}

	/**
	 * @return the chromosome
	 */
	public String getSequenceName() {
		return sequenceName;
	}

	/**
	 * @param sequenceName the chromosome to set
	 */
	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}

	@Override
	public boolean isNegativeStrand() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isPositiveStrand() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public int length() {
		// TODO Auto-generated method stub
		return last - first + 1;
	}

	/**
	 * @return the qualityScore
	 */
	public short getQualityScore() {
		return qualityScore;
	}

	/**
	 * @param qualityScore the qualityScore to set
	 */
	public void setQualityScore(short qualityScore) {
		this.qualityScore = qualityScore;
	}

	/**
	 * @param patternBeginning the patternBeginning to set
	 */
	public void setPatternBeginning(int patternBeginning) {
		this.patternBeginning = patternBeginning;
	}

	/**
	 * @param unitLength the unitLength to set
	 */
	public void setUnitLength(int unitLength) {
		this.unitLength = unitLength;
	}

	/**
	 * @param sumOfHeads the sumOfHeads to set
	 */
	public void setSumOfHeads(int sumOfHeads) {
		this.sumOfHeads = sumOfHeads;
	}

	/**
	 * @return the pattern
	 */
	public String getPattern() {
		return pattern;
	}

	/**
	 * @param pattern the pattern to set
	 */
	public void setPattern(String pattern) {
		this.pattern = pattern;
	}
	
	@Override
	public String toString() {
		return sequenceName + " " + first + " " + last + " " + unitLength + " " + numCopies + " " + length() + " " + pattern;
		
	}

}
