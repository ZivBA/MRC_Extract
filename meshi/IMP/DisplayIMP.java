package meshi.IMP;


import java.io.BufferedWriter;
import java.io.FileWriter;

public class DisplayIMP {

	private DomainList domainList;
	
	public DisplayIMP(DomainList domainList) {
		this.domainList = domainList;
	}
	
	// Main method
	public void writeToDisk(String projectName) {
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(projectName + "_display.xml"));
			bw.write("<Display>\n");
			for (DomainWithCrossLinks domain : domainList) {
				double red,green;
				switch (domain.proteinName().charAt(0)) {
		    	case 'A': 
		    		red=0;	
		    		green=0.7;
		    		break;
		    	case 'B': 
		    		red=0.2;	
		    		green=0.7;
		    		break;
		    	case 'G': 
		    		red=0.4;	
		    		green=0.7;
		    		break;
		    	case 'D': 
		    		red=0.6;	
		    		green=0.7;
		    		break;
		    	case 'E': 
		    		red=0;	
		    		green=0.3;
		    		break;
		    	case 'H': 
		    		red=0.2;	
		    		green=0.3;
		    		break;
		    	case 'Q': 
		    		red=0.4;	
		    		green=0.3;
		    		break;
		    	case 'Z': 
		    		red=0.6;	
		    		green=0.3;
		    		break;
		    	default:
		    		throw new RuntimeException("Invalid unit letter {A,B,G,D,H,Q,Z,W,E,I,J}");
		    	}
				bw.write("    <Protein id=\""+ domain.proteinName() + "_" + domain.domainName() +"\">\n");
				bw.write("        <Chain>\n");
				bw.write("            <Fragment id=\""+ domain.proteinName() + "_" + domain.domainName() +"_frag\">\n");
				bw.write("                <Color r=\""+red+"\" g=\""+green+"\" b=\""+(0.15+Double.valueOf(domain.domainName()).doubleValue()/2.5)+"\"/>\n");
				bw.write("            </Fragment>\n");
				bw.write("        </Chain>\n");
				bw.write("    </Protein>\n");
				for (CrossLinkedLysine xlK : domain.xlLysines() ) {
					bw.write("    <Protein id=\"" + domain.proteinName() + "_" + domain.domainName() + "_xl_" + xlK.resNum() + "\">\n");
					bw.write("        <Chain>\n");
					bw.write("            <Fragment id=\"" + domain.proteinName() + "_" + domain.domainName() + "_xl_" + xlK.resNum() + "_frag\">\n");
					bw.write("                <Color r=\""+1.0+"\" g=\""+1.0+"\" b=\""+1.0+"\"/>\n");
					bw.write("            </Fragment>\n");
					bw.write("        </Chain>\n");
					bw.write("    </Protein>\n");
				}
				for (CrossLinkedLysine xlK : domain.boundaries() ) {
					bw.write("    <Protein id=\"" + domain.proteinName() + "_" + domain.domainName() + "_bound_" + xlK.resNum() + "\">\n");
					bw.write("        <Chain>\n");
					bw.write("            <Fragment id=\"" + domain.proteinName() + "_" + domain.domainName() + "_bound_" + xlK.resNum() + "_frag\">\n");
					bw.write("                <Color r=\""+1.0+"\" g=\""+0.0+"\" b=\""+0.0+"\"/>\n");
					bw.write("            </Fragment>\n");
					bw.write("        </Chain>\n");
					bw.write("    </Protein>\n");
				}					
			}
			bw.write("</Display>\n");
			bw.close();
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}  		
	}

}
