package meshi.util;
import java.util.Iterator;
import java.util.StringTokenizer;

import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;
import meshi.util.string.StringList;

/**
 * Command File Parser.
 * <pre>
 * MESHI programs require a large number of user defined parameters (e.g. weights for the 
 * energy terms). The user provide these parameters through a single command file encoded in 
 * a very simple format. The CommandList class reads this file, parses it, stores the commands and 
 * provide them to the other classes.
 * 
 * The constructor of the class gets the command file name, and after the file has been successfully 
 * loaded many methods are available for the user to extract the parameters. 
 *
 * The rules for the command file format are:
 * empty lines - ignored
 * line starting with the "#" char - ignored, and used for commenting
 * other cases - the first word in the line serves as the key word 
 *
 *</pre>
 **/


public class CommandList extends MeshiList implements KeyWords { 
    /**
     * The file from which the comands were read. Usefull for debugging.
     **/
    private String comment;

    private boolean debug = true;
    
    /**
     * An empty commands List,
     **/
    public CommandList() {
	super(new IsCommand());
	this.comment = "commadList was not generated from a file";
    }

    public CommandList(CommandListGenerator generator) {
	this();
	this.comment  = generator.CommandListComment();
    }
    /**
     * Constract a CommandList object from a file.
     **/
    public CommandList(String commandsFileName) {
	this();
	this.comment = commandsFileName;
	StringList commands;
	MeshiLineReader commandsFile = new MeshiLineReader(commandsFileName);
	try {
	    commands = new StringList(commandsFile);
	} catch (Exception e) {throw new RuntimeException("A problem while reading commands file "+
							  commandsFileName+"\n"+e);}
	Iterator comanndsIter = commands.iterator();
	String line;
	while ((line = (String) comanndsIter.next()) != null) {
	    StringTokenizer tokenizer = new StringTokenizer(line);
	    if (tokenizer.hasMoreTokens() && (line.charAt(0) != '#'))
		add(new Command(line, comment));
	}
    }
    

    
    private void comment(String s) {
	comment = s;
    }

    public CommandList firstWordFilter(Key key) {
	return firstWordFilter(key.key);
    }

    /**
     * Generates a list including only the commands starting with a given keyword. 
     * Used to extract all the commands relevant to some module.
     **/  
    public CommandList firstWordFilter(String key) {
	CommandList out = new CommandList();
	out.comment(comment);
	Iterator iter = iterator();
	Command command;
	while ((command = (Command) iter.next()) != null)
	    if (command.firstWord().equals(key)) { 
		out.add(command);
		if (debug) command.printLine("====>");
	    }
	if (out.size() < 1) 
	    new RuntimeException("No command starts with "+key+" in\n"+
				 comment);
	return out;
    }
    /**
       Testing if some key is exists
     */
    public boolean keyExists(Key key){
	return keyExists(key.key);
    }

    /**
       Testing if some key is exists
     */
    public boolean keyExists(String key){
        for(Iterator i=iterator();i.hasNext();){
            Command c = (Command) i.next();
            if(c.firstWord().equals(key))
                return true;
        }
        return false;
    }

    /**
     * Returns the first command that start with the given keyword.
     **/
    public Command firstWord(Key key) {
	return  firstWord(key.key);
    }

    /**
     * Returns the first command that start with the given keyword.
     **/
    public Command firstWord(String key) {
	Iterator iter = iterator();
	Command command;

	while ((command = (Command) iter.next()) != null)
	    if (command.firstWord().equals(key)) {
		if (debug) command.printLine("=-=->");
		return command;
	    }
	throw new RuntimeException("No command starts with "+key+" in\n"+
				 comment);
    }

    /**
     * Returns the first command that whose second word is the given keyword. 
     **/
    public Command secondWord(Key key) {
	return secondWord(key.key);
    }

    /**
     * Returns the first command that whose second word is the given keyword. 
     **/
    public Command secondWord(String key) {
	Iterator iter = iterator();
	Command command;
	while ((command = (Command) iter.next()) != null)
	    if (command.secondWord().equals(key)) {
		if (debug) command.printLine("---->");
		return command;
	    }
	throw new RuntimeException("No command has "+key+" in its second position\n"+
				 comment);
    }

    /**
     * Reads a long string (protein sequence, secandary structure prediction ect.) that is divided into several lines.
     * Each line must start with the "keyword" followed by an integer. The reading is done with the 
     * word "end". Here is an example. In the command file you have:
     * myseq1 ACDEFGH
     * myseq2 IKLMN
     * myseq3 PQRSTVWY
     * myseq4 end
     * the method getSequence("myseq") will return "ACDEFGHIKLMPQRSTUVWY"
     **/
    public String getSequence(Key keyword) {
	return getSequence(keyword.key);
    }

    /**
     * Reads a long string (protein sequence, secandary structure prediction ect.) that is divided into several lines.
     * Each line must start with the "keyword" followed by an integer. The reading is done with the 
     * word "end". Here is an example. In the command file you have:
     * myseq1 ACDEFGH
     * myseq2 IKLMN
     * myseq3 PQRSTVWY
     * myseq4 end
     * the method getSequence("myseq") will return "ACDEFGHIKLMPQRSTUVWY"
     **/
    public String getSequence(String keyword) {
	Iterator iter = iterator();
	Command command;
	int i = 1;
	String key;
	String sequence = "";
	String fragment = "";
	while (! fragment.equals(END)) {
	    key = keyword+i;
	    command = firstWord(key);
	    fragment = command.secondWord();
	    if (! fragment.equals(END))
		sequence += fragment;
	    i++;
	}
	return sequence;
    }

    /**
     * Gets energy term weight getWeight.
     * If the following line occur in the command file: 
     * VDW weight 8.5 
     * then getWeight("VDW") returns 8.5
     **/
    public double getWeight(Key key) {
	return getWeight(key.key);
    }

    /**
     * Gets energy term weight getWeight.
     * If the following line occur in the command file: 
     * VDW weight 8.5 
     * then getWeight("VDW") returns 8.5
     **/
    public double getWeight(String key) {
	try {
	    return firstWordFilter(key).secondWord(WEIGHT).thirdWordDouble();
	} 
	catch (RuntimeException e) {
	    throw new RuntimeException("Failed to find weight for "+key+"\n"+e);
	}
   }

    private static class IsCommand implements Filter {
	public boolean accept(Object obj) { 
	    return (obj instanceof Command);
	}
    }
    
    public String comment() {return comment;}

    public void printCommandsLines(MeshiWriter writer) {
        for (int i = 0; i < size; i++) {
            writer.println(((Command)elementAt(i)).getLine());
        }
    }
}
    
