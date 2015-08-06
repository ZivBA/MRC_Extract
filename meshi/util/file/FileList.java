package meshi.util.file;
import java.io.File;
import java.util.Iterator;

import meshi.util.MeshiException;
import meshi.util.MeshiList;
import meshi.util.filters.Filter;
import meshi.util.string.StringList;

public class FileList extends MeshiList {
    //#2 constructors
    // generate a file list from either a text file with the names of the 
    // files or from a directory listing.
  public FileList() {
    super(new IsFile());
  }
  public FileList(File file) {
    this();
    if (file.isDirectory()) FileListFromDirectory(file);
    else FileListFromFile(file);
  }
  public FileList(String fileName) {
    this(new File(fileName));
  }
  public FileList(StringList fileNames) {
    this();
    addFormStringList(fileNames);
  }
  //---------------------------------------------------------------------------
    // #3a constructor supplaments
    private void FileListFromFile(File file) {
	StringList fileNames = 
	    (new StringList(new MeshiLineReader(file))).flatten();
	addFormStringList(fileNames);
    }
    private void FileListFromDirectory(File file) {
	throw new MeshiException("FileList.fileListFromDirectory(File file) error:"+
			      "option not yet supported\n"+
			      "file: "+file+"\n");
    }
    public void addFormStringList(StringList fileNames) {
	Iterator iter = fileNames.iterator();
	String fileName;
	while ((fileName = (String) iter.next()) != null)
	    add(new File(fileName));
    }
     static class IsFile implements Filter {
	public boolean accept(Object obj) {
		return (obj instanceof File);
	}
    }
   
    //-------------------------------------------------------------------------
    public boolean sortable() {return false;}
}

    
