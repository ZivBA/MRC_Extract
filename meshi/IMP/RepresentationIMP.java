package meshi.IMP;

import java.io.BufferedWriter;
import java.io.FileWriter;

public class RepresentationIMP {

	private DomainList domainList;
	
	public RepresentationIMP(DomainList domainList) {
		this.domainList = domainList;
	}
	
	// Main method
	public void writeToDisk(String projectName) {
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(projectName + "_representation.xml"));
			bw.write("<Representation>\n");
			for (DomainWithCrossLinks domain : domainList) {
				bw.write("    <Protein id=\""+ domain.proteinName() + "_" + domain.domainName() +"\">\n");
				bw.write("        <Chain> id=\""+ domain.proteinName() + "_" + domain.domainName() +"_chain\">\n");
				bw.write("            <Fragment id=\""+ domain.proteinName() + "_" + domain.domainName() +"_frag\">\n");
				bw.write("                <GeometricShapeRep total_residue=\"" + domain.numberOfResidues() + "\">\n");
				bw.write("                    <Sphere><InitialPosition x=\""+domain.x()+"\" y=\""+domain.y()+"\" z=\""+domain.z()+"\" optimize=\"1\"/></Sphere>\n");
				bw.write("                </GeometricShapeRep>\n");
				bw.write("            </Fragment>\n");
				bw.write("        </Chain>\n");
//				bw.write("    </Protein>\n");
				for (CrossLinkedLysine xlK : domain.xlLysines() ) {
//					bw.write("    <Protein id=\""+ domain.proteinName() + "_" + domain.domainName() + "_xl_" + xlK.resNum() + "\">\n");
					bw.write("        <Chain id=\""+ domain.proteinName() + "_" + domain.domainName() + "_xl_" + xlK.resNum() + "_chain\">\n");
					bw.write("            <Fragment id=\""+ domain.proteinName() + "_" + domain.domainName() + "_xl_" + xlK.resNum() + "_frag\">\n");
					bw.write("                <GeometricShapeRep total_residue=\"90\">\n");
					bw.write("                    <Sphere><InitialPosition x=\""+xlK.x()+"\" y=\""+xlK.y()+"\" z=\""+xlK.z()+"\" optimize=\"1\"/></Sphere>\n");
					bw.write("                </GeometricShapeRep>\n");
					bw.write("            </Fragment>\n");
					bw.write("        </Chain>\n");
//					bw.write("    </Protein>\n");
				}
				for (CrossLinkedLysine xlK : domain.boundaries() ) {
//					bw.write("    <Protein id=\"" + domain.proteinName() + "_" + domain.domainName() + "_bound_" + xlK.resNum() + "\">\n");
					bw.write("        <Chain id=\"" + domain.proteinName() + "_" + domain.domainName() + "_bound_" + xlK.resNum() + "_chain\">\n");
					bw.write("            <Fragment id=\"" + domain.proteinName() + "_" + domain.domainName() + "_bound_" + xlK.resNum() + "_frag\">\n");
					bw.write("                <GeometricShapeRep total_residue=\"30\">\n");
					bw.write("                    <Sphere><InitialPosition x=\""+xlK.x()+"\" y=\""+xlK.y()+"\" z=\""+xlK.z()+"\" optimize=\"1\"/></Sphere>\n");
					bw.write("                </GeometricShapeRep>\n");
					bw.write("            </Fragment>\n");
					bw.write("        </Chain>\n");
//					bw.write("    </Protein>\n");					
				}
				bw.write("    </Protein>\n");
			}
			bw.write("</Representation>\n");
			bw.close();
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}  		
	}

}
