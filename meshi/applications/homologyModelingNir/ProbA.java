package meshi.applications.homologyModelingNir;

import java.util.StringTokenizer;

import meshi.util.file.File2StringArray;

public class ProbA {
	
	private String[] suboptimalTops;
	private String[] suboptimalBottoms;
	private double[] scores;
	private String top;
	private String bottom;
	
	public ProbA(String probAfilename) {
		StringTokenizer st;
		
		// Reading the sequences
		String[] content = File2StringArray.f2a(probAfilename);
		int lineNum = 0;
		while (!content[lineNum].startsWith("#upper"))
			lineNum++;
		lineNum+=2;
		st = new StringTokenizer(content[lineNum]);
		st.nextToken(); // Reading the first #
		top = st.nextToken();
		System.out.println("Upper Sequence:\n" + top);
		lineNum+=3;
		st = new StringTokenizer(content[lineNum]);
		st.nextToken(); // Reading the first #
		bottom = st.nextToken();
		System.out.println("Lower Sequence:\n" + bottom);
		lineNum+=2;
		
		// Reagong the alignments
		int alignmentCounter = 0;
		suboptimalBottoms = new String[content.length-lineNum];
		suboptimalTops = new String[content.length-lineNum];
		scores = new double[content.length-lineNum];
		for (;lineNum<content.length ; lineNum++) {
			st = new StringTokenizer(content[lineNum]);
			String alignment = st.nextToken();
			scores[alignmentCounter] = (new Double(st.nextToken())).doubleValue(); 
			suboptimalBottoms[alignmentCounter] = "";
			suboptimalTops[alignmentCounter] = "";
			int topCounter = 0;
			int bottomCounter = 0;
			for (int inLineCounter=0 ; inLineCounter<alignment.length() ; inLineCounter++) {
				if (alignment.charAt(inLineCounter)=='.') {
					suboptimalBottoms[alignmentCounter] += '-';
					suboptimalTops[alignmentCounter] += top.charAt(topCounter);
					topCounter++;
				}
				else if (alignment.charAt(inLineCounter)=='^') {
					suboptimalBottoms[alignmentCounter] += bottom.charAt(bottomCounter);
					suboptimalTops[alignmentCounter] += '-';
					bottomCounter++;
				}
				else {
					suboptimalTops[alignmentCounter] += top.charAt(topCounter);
					suboptimalBottoms[alignmentCounter] += bottom.charAt(bottomCounter);
					topCounter++;
					bottomCounter++;					
				}				
			}
			alignmentCounter++;
		}		
	}

	public int alignmentNum() {
		return scores.length;
	}
	
	public String suboptimalTops(int i) {
		return suboptimalTops[i];
	}

	public String suboptimalBottoms(int i) {
		return suboptimalBottoms[i];
	}

	public double scores(int i) {
		return scores[i];
	}

	public String getTop() {
		return top;
	}

	public String getBottom() {
		return bottom;
	}
	
	public String toString() {
		String out = "";
		for (int c=0 ; c<alignmentNum() ; c++) {
			out += (scores(c) + "\n");
			out += (suboptimalTops(c) + "\n");
			out += (suboptimalBottoms(c) + "\n");
			out += "\n";
		}
		return out;		
	}

}
