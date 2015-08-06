package meshi.IMP;

import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.util.crossLinking.Crosslink;
import meshi.util.crossLinking.CrosslinkVector;

public class RestraintIMP {

	private DomainList domainList;
	
	public RestraintIMP(DomainList domainList) {
		this.domainList = domainList;
	}
	
	// Main method
	public void writeToDisk(String projectName , CrosslinkVector xlVecRealXLs, CrosslinkVector xlVecBoundaries) {
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(projectName + "_restraints.xml"));
			bw.write("<RestraintSet>\n");
			// Excluded Volume
			// ---------------
			bw.write("    <ExcludedVolume>\n");
			bw.write("        <Restraint>\n");
			for (DomainWithCrossLinks domain : domainList) {
				bw.write("            <Particle id=\""+ domain.proteinName() + "_" + domain.domainName() +"_frag\"/>\n");
			}
			bw.write("        </Restraint>\n");			
			bw.write("    </ExcludedVolume>\n");
			// Rigid Body
			// ----------
			bw.write("    <RigidBody>\n");
			for (DomainWithCrossLinks domain : domainList) {
				bw.write("        <Restraint>\n");
				bw.write("            <Particle id=\""+ domain.proteinName() + "_" + domain.domainName() +"\"/>\n");
//				for (CrossLinkedLysine xlK : domain.xlLysines() ) {
//					bw.write("            <Particle id=\""+ domain.proteinName() + "_" + domain.domainName() +"_xl_" + xlK.resNum() + "\"/>\n");					
//				}
//				for (CrossLinkedLysine xlK : domain.boundaries() ) {
//					bw.write("            <Particle id=\""+ domain.proteinName() + "_" + domain.domainName() +"_bound_" + xlK.resNum() + "\"/>\n");					
//				}				
				bw.write("        </Restraint>\n");			
			}
			bw.write("    </RigidBody>\n");
			// Y2H
			// ---
			bw.write("    <Y2H>\n");
			for (Crosslink xl : xlVecRealXLs) {
				Domain dom1 = domainList.findDomain(xl.protName1(), xl.absPos1());
				Domain dom2 = domainList.findDomain(xl.protName2(), xl.absPos2());
				if ((dom1!=null) && (dom2!=null)) {
					bw.write("        <Restraint>\n");
					bw.write("            <Particle id=\""+ xl.protName1() + "_" + dom1.domainName() + "_xl_" + xl.absPos1() + "_frag\"/>\n");
					bw.write("            <Particle id=\""+ xl.protName2() + "_" + dom2.domainName() + "_xl_" + xl.absPos2() + "_frag\"/>\n");				
					bw.write("        </Restraint>\n");							
				}
			}			
			for (Crosslink xl : xlVecBoundaries) {
				Domain dom1 = domainList.findDomain(xl.protName1(), xl.absPos1());
				Domain dom2 = domainList.findDomain(xl.protName2(), xl.absPos2());
				if ((dom1!=null) && (dom2!=null)) {
					bw.write("        <Restraint>\n");
					bw.write("            <Particle id=\""+ xl.protName1() + "_" + dom1.domainName() + "_bound_" + xl.absPos1() + "_frag\"/>\n");
					bw.write("            <Particle id=\""+ xl.protName2() + "_" + dom2.domainName() + "_bound_" + xl.absPos2() + "_frag\"/>\n");				
					bw.write("        </Restraint>\n");							
				}
			}			
			bw.write("    </Y2H>\n");			
			bw.write("</RestraintSet>\n");				
			bw.close();
		}
		catch(Exception e) {
			System.out.print(e);
			throw new RuntimeException(e.getMessage());
		}  		
	}

}
