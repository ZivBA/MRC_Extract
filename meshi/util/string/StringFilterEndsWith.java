package meshi.util.string;
public  class StringFilterEndsWith extends StringFilter {
    public  StringFilterEndsWith(String key) {super(key);}
    public  StringFilterEndsWith(StringList keys) {super(keys);}
    public boolean accept(String string,String key) {
	return string.endsWith(key);
    }
}
	
    
