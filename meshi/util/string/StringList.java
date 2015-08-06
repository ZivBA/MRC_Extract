package meshi.util.string;
import java.io.File;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.util.Iterator;
import java.util.StringTokenizer;

import meshi.util.MeshiException;
import meshi.util.MeshiList;
import meshi.util.file.MeshiLineReader;
import meshi.util.filters.Filter;


// constructors
public class StringList extends MeshiList {
    public StringList() {
	super(new IsString());
	if (filter() == null) throw new RuntimeException("filter is null");
    }
    public StringList(MeshiListIterator iterator) throws Exception{
        super(iterator);
    }
    public StringList(String string) {
	this();
	add(string);
    }
    public StringList(StringTokenizer st) {
	    this();
	    while (st.hasMoreTokens()) add(st.nextToken());
    }
    public StringList(String parseMe, StringList separators) {
	this();
	String string;
	Iterator parseIter = 
	    (StringParser.bySeparators(parseMe,
					 separators)).iterator();
	while ((string = (String) parseIter.next()) != null)
	    add(string);
    }
     public StringList(String parseMe, String separator) {
	this(parseMe,new StringList(separator));
    }
    public StringList(StreamTokenizer ST) {
	this();
	ST.wordChars('_','_');
	ST.wordChars('/','/');
	ST.wordChars('1','1');
	try
	    {
		while(ST.nextToken() != StreamTokenizer.TT_EOF)
		    {
			if (ST.sval != null) add(ST.sval);
			else add((new Double(ST.nval)).toString());
		    }
	    }
	catch(IOException e)
                { 
		    throw new MeshiException("StringList(StreamTokenizer ST) Error:\n"+
				       e.getMessage());
		}
    }
    public StringList(MeshiLineReader MLR) {
	this();
	setComment(MLR.path());
	String string;
	try{
	    while((string = MLR.readLine()) != null)
		add(string);
	    MLR.close();
	}
	catch(Exception e){ 
	    throw new MeshiException("StringList(MeshiLineReader MLR) Error:\n"+
				     e.getMessage());
	}
    }

    public StringList(File file) {
	this(new MeshiLineReader(file));
    }

    public StringList(MeshiLineReader MLR, Filter filter) {
	this();
	setComment(MLR.path());
	String string;
	try
	    {
		while((string = MLR.readLine()) != null)
		    if (filter.accept(string)) add(string);
	    }
	catch(Exception e)
                { 
		    throw new MeshiException("StringList(MeshiLineReader MLR) Error:\n"+
				       e.getMessage());
		}
    }
	    
			
    public String stringAt(int index) { return (String) elementAt(index);}
    public String lastString() { return stringAt(size() - 1);}
    public StringList filter(Filter filter) {
	return (StringList) filter(filter,new StringList());
    }
 
//     public StringList filter(StringFilter filter) {
// 	StringList newList = new StringList();
// 	String string;
// 	Iterator listIterator = iterator();
// 	while ((string = (String) listIterator.next())!= null)
// 	    if (filter.accept(string)) newList.add(string);
// 	return newList;
//     }
                  
    public String startsWith(String key) {
	String string;
	Iterator listIterator = iterator();
	while ((string = (String) listIterator.next())!= null) {
	    if (string.startsWith(key)) return string;
	}
	return null;
    }

    public StringList filterStartsWith(String key) {
	StringFilterStartsWith filter = 
	    new StringFilterStartsWith(key);
	return filter(filter);
    }
    public StringList filterStartsWith(StringList keys) {
	StringFilterStartsWith filter = 
	    new StringFilterStartsWith(keys);
	return filter(filter);
    }
    public StringList filterEndsWith(String key) {
	StringFilterEndsWith filter = 
	    new StringFilterEndsWith(key);
	return filter(filter);
    }
    public StringList filterEndsWith(StringList keys) {
	StringFilterEndsWith filter = 
	    new StringFilterEndsWith(keys);
	return filter(filter);
    }
    public StringList filterGrep(String key) {
	StringFilterGrep filter = 
	    new StringFilterGrep(key);
	return filter(filter);
    }
    public StringList filterGrep(StringList keys) {
	StringFilterGrep filter = 
	    new StringFilterGrep(keys);
	return filter(filter);
    }
    public StringList stringParseAt(int index,StringList separators) {
	return StringParser.bySeparators(stringAt(index),separators);
    }
    public StringList stringParseAt(int index,String separator) {
	return StringParser.bySeparator(stringAt(index),separator);
    }
    public StringList stringParseAt(int index) {
	return StringParser.standard(stringAt(index));
    }
    public StringList flatten() {
	StringList newList = new StringList(),tempList;
	Iterator SI = iterator();
	Iterator SI1;
	String string;
	while ((tempList = StringParser.standard((String) SI.next())) != null)
	    {
		SI1 = tempList.iterator();
		while ((string = (String) SI1.next()) != null)
		    newList.add(string);
	    }
	return newList;
    }
    
	
    public static StringList standardSeparators() {
	StringList newList = new StringList();
	newList.add(" ");
	newList.add(",");
	newList.add(";");
	newList.add("\n");
	newList.add("\t");
	return newList;
    }
    static class IsString implements Filter {
	public boolean accept(Object obj) {
		return (obj instanceof String);
	}
    }
    public boolean sortable() {return true;}
    public StringList filterStartsWithKeyAndCut(StringList keysList) {
        StringList ans = new StringList();
        Iterator keys = keysList.iterator();
        String key;
        boolean foundKey = false;
        Iterator lines = this.iterator();
        String line;
        while ((line = (String)lines.next())!=null){
            while ((!foundKey)&&((key = (String)keys.next()) != null)){
              if (line.startsWith(key)){
                 foundKey = true;
                 ans.add(line.substring(key.length()));
              }
            }
            foundKey = false;
            keys = keysList.iterator();
        }
        return ans;
    }

}	
		
	
