package stats;

import java.io.File;

import tools.Log;

public class StatsCalc {

	/** Hyper-Parallelized Composite of Multiple Signals (CMS) Java implementation GWAS and Local study version 1.0
	 * This program calculated the stats that are used in CMS analysis
	 * @author Hayden Smith
	 * 
	 * @param args0		out dir
	 * @param args1		chr
	 * @param args2		window number
	 * @param args3		threading
	 */
	public static void main(String[] args) {

		Log log = new Log();
		log.addLine("\t\t\t*****Starting CMS Stat Calculation*****\n");
		
		File out_dir = new File(args[0]);
		if(!out_dir.exists() && !out_dir.isDirectory()) {
			//TODO: throw some error
		}
		
		int chr = Integer.parseInt(args[1]);
		int win_num = Integer.parseInt(args[2]);
		boolean threading = Boolean.parseBoolean(args[3]);
		
		StatsCalc sc = new StatsCalc(out_dir, chr, win_num, threading, log);
		sc.runStats();
	}
	
	private boolean threading;
	private int chr;
	private int win_num;
	
	private File out_dir;
	
	private Log log;
	
	public StatsCalc(File out_dir, int chr, int win_num, boolean threading, Log log) {
		
		this.threading = threading;
		this.chr = chr;
		this.win_num = win_num;
		this.out_dir = out_dir;
		this.log = log;
	}
	
	public void runStats() {
		
		//create a stats_file dir inside of out_dir
		//load envi variables
		//create all the calculators
		//run the tests
		//write the window stats output
	}
	
	

}
