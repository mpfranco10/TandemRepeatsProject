package uniandes.algorithms.tr;


import java.io.FileWriter;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;

import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;

public class TandemRepeatsExample {

	public static void main(String[] args) throws Exception {
		
		
		String fastaFilename = "./data/S288C_20150113.fa";
		// Assemble sequence
		

		FastaSequencesHandler handler = new FastaSequencesHandler();
		handler.setSequenceType(StringBuilder.class);
		QualifiedSequenceList sequences = handler.loadSequences(fastaFilename);
		if (sequences.size() == 0)
			throw new Exception("No sequences found in file: " + fastaFilename);
		Files.deleteIfExists(Paths.get("data\\outTR.txt")); 
		PrintWriter pw = new PrintWriter(new FileWriter("data\\outTR.txt", true));
		String line = "sequenceName Start End Period NumReps TotalSize Pattern Sequence IdealSeq Alignment AlignmentScore";
		pw.println(line);
		pw.close();
		
		long timeT = System.currentTimeMillis();
		for (int i = 0; i < sequences.size() ; i++) {
			
			QualifiedSequence seq = sequences.get(i);
			String sequence = seq.getCharacters().toString().toUpperCase();
			int seqLength = sequence.length();
			System.out.println("Length of the sequence read: " + seqLength);
			
			double mProb = 0.8;
			TRFCandidateSelector trfc = new TRFCandidateSelector(sequence, seq.getName(), mProb, 0.1, 25, 30, 2, 4);
			
			
				ArrayList<TandemRepeat> candidates = trfc.getAllCandidates();
				
				long time = System.currentTimeMillis();
				trfc.allignCandidates(candidates);
				time = System.currentTimeMillis()-time;
				time = time / 1000;
				System.out.println("Time performing allignments: "+time + " secs");
			
			
		}
		timeT = System.currentTimeMillis()-timeT;
		timeT = timeT / 1000;
		System.out.println("Total time: "+timeT + " secs");
	}

}
