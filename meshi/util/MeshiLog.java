package meshi.util;
import java.util.Formatter;

import meshi.util.string.StringList;


public class MeshiLog extends StringList {
    public final int LINE_LENGTH = 240; 
    public final String KEY = "MESHILOG"; 
    public final int FIELD_LENGTH = 10;
    public final int FIELDS_LENGTH = 6;

    public MeshiLog() { 
	super(); 
    } 

    public void add(String log, String key) {
	String formatString = KEY+" "+key+" %s";
	int effectiveLineLength = LINE_LENGTH - KEY.length() - key.length()-7;
	String line;
	while (log.length() > 0) { 
	    if (log.length() < effectiveLineLength) {
		Formatter frmt = new Formatter();
		frmt.format(formatString,log);
		line = frmt.out().toString();
		log = ""; 
	    } 
	    else { 
		int end = effectiveLineLength;
		while ((log.charAt(end) != ' ') & (end > 0))end--; 
		if (end == 0) throw new RuntimeException("bad formated log "+log);
		Formatter frmt = new Formatter();
		frmt.format(formatString,log.substring(0,end+1));
		line = frmt.out().toString();
		log = log.substring(end+1);
	    } 
	    add(line);
	} 
    }



    public String field(String s) {
	    return field(s,FIELD_LENGTH);
    }
    public String field(double d) {
	    return field(d,FIELD_LENGTH);
    }
    public String fields(String s) {
	    return field(s,FIELDS_LENGTH);
    }
    public String fields(double d) {
            return field(d,FIELDS_LENGTH);

    }

    public String field(String s,int length) {
	String out;
	if (s.length() <= length-1) {
	    out = s;
	    for (int i = 0; i <  (length-s.length());i++)
		out +=" ";
	}
	else out = s.substring(0,length-1)+" ";
	return out;
    }

    public String field(double d, int length) {
	String fmt;
	if (length == FIELD_LENGTH)  
	      fmt = "%-"+(length-1)+".2f";
	else 
   	      fmt = "%-"+(length-1)+".3f";
	String dString = (new Formatter()).format(fmt, d).toString();
	if ((dString.indexOf(".") == -1) &
	    (dString.length() >= length)) {
	    String dummy = "";
	    for (int i = 0; i < (length-1); i++) dummy += "*";
	    dummy += " ";
	    return dummy;
	}
	else return dString+" ";
    }	    
}
	



			
