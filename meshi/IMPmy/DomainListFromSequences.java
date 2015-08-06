package meshi.IMPmy;

import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.util.crossLinking.Crosslink;
import meshi.util.crossLinking.CrosslinkVector;
import meshi.util.crossLinking.FastaSeq;
import meshi.util.crossLinking.MySequence;
import meshi.util.crossLinking.MySequenceList;

public class DomainListFromSequences extends DomainList {

	public final int PARSING_SEQ_LENGTH = 10;	
	
	public  DomainListFromSequences(String fastaFileName) {
		MySequenceList seqList = new MySequenceList(fastaFileName);
		for (MySequence seq : seqList) {
			int domainCounter = 0;
			for (int startRes=0; startRes<seq.seq().length() ; startRes+=PARSING_SEQ_LENGTH) {
				int lastRes = Math.min(startRes+PARSING_SEQ_LENGTH-1, seq.seq().length());
				int[][] residues = {{startRes , lastRes}};
				Domain dom = new Domain(seq.title(), ""+domainCounter, residues, false);
				dom.setCenter(new AnchorPosition(seq.title(), ""+domainCounter, (startRes + lastRes)/2, 0.0, 0.0, 0.0));
				add(dom);
				domainCounter++;
			}
			System.out.println("Protein: " + seq.title() + " was parsed into " + domainCounter + " domains.");
		}
	}
	
	@Override
	public void setXLvecs(CrosslinkVector xlVec, double xlWeight) {
		for (Crosslink xl : xlVec) {
			Domain dom1 = findDomain(xl.protName1(), xl.absPos1());			
			Domain dom2 = findDomain(xl.protName2(), xl.absPos2());
			if (dom1==null) {
				throw new RuntimeException("\n\nCould not find domain for residue " + xl.absPos1() + " in protein " + xl.protName1() + "\n\n");
			}
			if (dom2==null) {
				throw new RuntimeException("\n\nCould not find domain for residue " + xl.absPos2() + " in protein " + xl.protName2() + "\n\n");
			}
			if (dom1!=dom2) {
				AnchorPosition pos1 = dom1.center();
				AnchorPosition pos2 = dom2.center();
				disConstList().add(new DistanceConstraint(pos1, pos2, DistanceConstraintType.CROSS_LINK,xlWeight));
			}
		}
	}

	@Override
	public void setBoundaries(double boundaryWeight) {
		for (int c=1 ; c<size() ; c++) {
			if (get(c).proteinName().equals(get(c-1).proteinName())) {
				disConstList().add(new DistanceConstraint(get(c).center(), get(c-1).center(), DistanceConstraintType.CONNECTIVITY, boundaryWeight));
			}
		}
	}

	@Override
	public void rigidify(double rigidityWeight) {
		throw new RuntimeException("Do not invoke Rigidify in this type of representation.");
	}

	@Override
	public void addEV(double evWeight) {
		for (int c1=0 ; c1<size() ; c1++) {
			for (int c2=c1+1 ; c2<size() ; c2++) {
				disConstList().add(new DistanceConstraint(get(c1).center(), get(c2).center(), DistanceConstraintType.EXCLUDED_VOLUME,evWeight));
			}
		}
	}

	@Override
	public void report(int reportNumber) {
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter("report_"+reportNumber+".txt"));			
			for (Domain dom : this) {
				bw.write(dom.proteinName() + " " + dom.domainName() + " " + dom.center().x() + " " + dom.center().y() + " " + dom.center().z() + "\n");
			}
			
			bw.close();
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}  				
	}
	
	
	public void seperateProteins(double boundingBox) {
		double x = boundingBox*(Math.random()-0.5);
		double y = boundingBox*(Math.random()-0.5);
		double z = boundingBox*(Math.random()-0.5);
		double counter = 0.0;
		for (int c=0 ; c<size()-1 ; c++) {
			counter += 2.0;
			get(c).moveCenterTo(x+counter, y+Math.random(), z+Math.random());
			if (!get(c).proteinName().equals(get(c+1).proteinName())) {
				x = boundingBox*(Math.random()-0.5);
				y = boundingBox*(Math.random()-0.5);
				z = boundingBox*(Math.random()-0.5);
				counter = 0.0;			
			}
		}
		counter += 2.0;
		lastElement().moveCenterTo(x+counter, y+Math.random(), z+Math.random());
	}

}
