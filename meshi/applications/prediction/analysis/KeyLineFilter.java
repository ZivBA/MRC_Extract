package meshi.applications.prediction.analysis;
import java.util.StringTokenizer;

public class KeyLineFilter extends  MeshiLogLineFilter {
    public final String key;
    public KeyLineFilter(String key) {
	this.key = key;
    }

    public boolean accept(Object obj) {
	StringTokenizer tokenizer = acceptLine((String) obj);
	if (tokenizer == null) return false;

	if (!tokenizer.hasMoreTokens()) return false;
	String token = tokenizer.nextToken();
	if (token.equals(key)) return true;
	return false;
	}
}
