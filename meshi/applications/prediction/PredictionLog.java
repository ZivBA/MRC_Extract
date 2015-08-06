package meshi.applications.prediction;
import java.util.Iterator;

import meshi.energy.AbstractEnergy;
import meshi.energy.TotalEnergy;
import meshi.energy.solvate.SolvateEnergy;
import meshi.util.MeshiLog;


public class PredictionLog extends MeshiLog {
    private TotalEnergy energy = null; 
    private double rms0 = -1, rms1 = -1, rms2 = -1; 
    private double gdt0 = -1, gdt1 = -1, gdt2 = -1; 
    public PredictionLog() { 
	super(); 
    } 

    public void setEnergy(TotalEnergy energy) {
	this.energy = energy;
    }

    public void setRms0(double rms0) { this.rms0 = rms0;  }

    public void setRms1(double rms1) {this.rms1 = rms1;}

    public void setRms2(double rms2) {this.rms2 = rms2;}

    public void setGdt0(double gdt0) { this.gdt0 = gdt0;  }

    public void setGdt1(double gdt1) {this.gdt1 = gdt1;}

    public void setGdt2(double gdt2) {this.gdt2 = gdt2;}

    public void comments() {
	String out = "";
	add( "rms0 - RMS of the initial model before loops are added and before any refinement.", "COMMENT");
	add( "rms1 - total RMS", "COMMENT");
	add( "rms2 - RMS of residues directly extracted from the template", "COMMENT");
	add( "gdt0 - GDT of the initial model before loops are added and before any refinement.", "COMMENT");
	add( "gdt1 - total GDT", "COMMENT");
	add( "gdt2 - GDT of residues directly extracted from the template", "COMMENT");
	add( "cold1 - total energy of the cold atoms (avg+std)", "COMMENT");
	add( "cold2 - total energy of the cold atoms (avg+2std)", "COMMENT");
	add( "Tn - titles of set #n", "COMMENT");
	add( "Vn - Values of set #n", "COMMENT");
	add( "", "COMMENT");
	add( "", "COMMENT");
	add( "", "COMMENT");
                                                                                                                                            
    }
   
    public void summary1() {
	    summary0("T1","V1");
    }
    public void summary3() {
	    energy.off();
            AbstractEnergy solvate = energy.getEnergyTerm(new SolvateEnergy());
	    solvate.on();
	    summary0("T3","V3");
    }
    public void summary0(String t, String v) {
	String titles = "";
	String values = "";
	
	titles += field("rms1");
	values += field(rms1);
	titles += field("rms2");
	values += field(rms2);
	titles += field("gdt1");
	values += field(gdt1);
	titles += field("gdt2");
	values += field(gdt2);
	
	titles += field("TotalEnergy");
	try {
		energy.update();
		energy.evaluate();
		energy.evaluateAtoms();
	}
	catch (Exception ex) {throw new RuntimeException("Exception in log\n"+ex);}
	values += field(energy.getLastEnergy());
	titles += field("avgEnergy");
	values += field(energy.avgEnergy());
	titles += field("cold1");
	values += field(energy.filteredEnergy(1));
	titles += field("cold2");
	values += field(energy.filteredEnergy(2));
	titles += field("avgCold1");
	values += field(energy.avgFilteredEnergy(1));
	titles += field("avgCold2");
	values += field(energy.avgFilteredEnergy(2));
	
	add(titles,t);
	add(values,v);
    }	
    
    public void summary2() {
	String titles = "";
	String values = "";
	
	titles += field("rms1");
	values += field(rms1);
	titles += field("rms2");
	values += field(rms2);
	titles += field("gdt1");
	values += field(gdt1);
	titles += field("gdt2");
	values += field(gdt2);
	titles += field("TotalEnergy");
	values += field(energy.getLastEnergy());
	titles += field("avgEnergy");
	values += field(energy.avgEnergy());
	
	Double value;
	Iterator energyValues = energy.energyValues().iterator();;
	Iterator terms = energy.energyTerms().iterator();
	while ((value = (Double)energyValues.next()) != null) {
            String name = ((AbstractEnergy) terms.next()).comment();
            titles += field(name);
	    values += field(value.doubleValue());
        }
	add(titles,"T2");
	add(values,"V2");
    }
    public void summary5(int numberOfAtoms) {
        try {
		energy.on();
                energy.update();
                energy.evaluate();
                energy.evaluateAtoms();
            }
           catch (Exception ex) {throw new RuntimeException("Exception in log\n"+ex);}

	String titles = "";
	String values = "";
	
	titles += fields("rms0");
	values += fields(rms0);
	titles += fields("rms1");
	values += fields(rms1);
	titles += fields("rms2");
	values += fields(rms2);
	titles += fields("gdt0");
	values += fields(gdt1);
	titles += fields("gdt0");
	values += fields(gdt1);
	titles += fields("gdt2");
	values += fields(gdt2);
        titles += field("TotalEnergy");
        values += field(energy.getLastEnergy());
        titles += fields("avgEnergy");
        values += fields(energy.avgEnergy());

	String[] keys = {"Solvate","2Torsion","HB-PairsE","propen"};
	Double value;
	for (int i = 0 ; i < keys.length; i++) {
		Iterator energyValues = energy.energyValues().iterator();
		Iterator terms = energy.energyTerms().iterator();

	    String key = keys[i];
	    while ((value = (Double)energyValues.next()) != null) {
		String name = ((AbstractEnergy) terms.next()).comment();
		if (name.startsWith(key)) {
		    titles += field(name);
		    values += field(value.doubleValue());
		    titles += fields("Avg"+name);
		    values += fields(value.doubleValue()/numberOfAtoms);
		}
	    }
	}
	add(titles,"T5");
	add(values,"V5");
    }
    
    public void summary4() {
	    String titles = "";
	    String values = "";

	    titles += field("rms0");
	    titles += field("rms1");
	    titles += field("rms2");
	    values += field(rms0);
	    values += field(rms1);
	    values += field(rms2);
	    titles += field("gdt0");
	    titles += field("gdt1");
	    titles += field("gdt2");
	    values += field(gdt0);
	    values += field(gdt1);
	    values += field(gdt2);

	    add(titles,"T4");
	    add(values,"V4");
    }
}
	



			
