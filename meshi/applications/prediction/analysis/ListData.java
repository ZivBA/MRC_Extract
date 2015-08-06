package meshi.applications.prediction.analysis;
import java.util.Formatter;

public class ListData {
    public final double avg, std;
    public final int number;
    public final String comment;
    public ListData(double avg, double std, int number, String comment) {
	this.avg = avg;
	this.std = std;
	this.number = number;
	this.comment = comment;
    }
    public String toString() {
	StringBuffer res = new StringBuffer();
	Formatter f = new Formatter(res); 
	f.format("%8.3f(%-6.3f) %-5d %s",avg,std, number, comment);
	f.flush();
	f.close();
	return res.toString();
    }
}
