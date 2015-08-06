package meshi.applications.prediction.analysis;
import java.util.StringTokenizer;

import meshi.util.KeyWords;
import meshi.util.filters.Filter;

public class MeshiLogLineFilter implements Filter, KeyWords {
    public static StringTokenizer acceptLine(String line) {
	StringTokenizer tokenizer = new StringTokenizer(line);
	if (! tokenizer.hasMoreTokens()) return null;
	String token = tokenizer.nextToken();
	if (! token.equals(MESHILOG_KEY.key)) return null;
	return  tokenizer;
    }
    
    public boolean accept(Object obj) {
	String line = (String) obj;
	StringTokenizer tokenizer = acceptLine(line);
	if (tokenizer == null) return false;
	return true;
	}
}
