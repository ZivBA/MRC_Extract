package meshi.util.string;
public  class StringFilterStartsWith extends StringFilter {
    public  StringFilterStartsWith(String key) {super(key);}
    public  StringFilterStartsWith(StringList keys) {super(keys);}
    public boolean accept(String string,String key) {
	return string.startsWith(key);
    }
}
	
    
