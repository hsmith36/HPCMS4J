package analysis;

import java.io.File;
import java.util.LinkedList;
import java.util.List;

import errors.IllegalInputException;
import tools.*;

public class Combiner {

	/**
	 * Gets all the data from all the stats output files
	 * If directly called from main it will print the data (any combination of tests) into one big file
	 * If called from Analyzer it will return List<WindowStats> for further analysis
	 * 
	 * @param args0		out dir
	 * @param args1		stats string ("i:h:x:d:f" or any combination of that)
	 * @param args2		window range ("x-y")
	 */
	public static void main(String[] args) {
		
		Log log = new Log(Log.type.combine);
		
		try {
			 checkArgs(args, log);
			 
			 Combiner c = new Combiner(args[0], log);
			 c.combineStats(args[1], args[2]);
			 c.writeStats();
			 
			 log.addLine("Finished!");
			 
		} catch (Exception e) {
			
			System.out.println("CMS Died Prematurely." 
					+ " Check log output for troubleshooting.");
			
			log.addLine("\n\nCMS Died Prematurely. Error in computation.");
			
			e.printStackTrace();
		}
	}
	
	private static void checkArgs(String[] args, Log log) 
			throws IllegalInputException {
		
		if(args.length != 3) {
			String msg = "Error: Parameter length incorrect";
			throw new IllegalInputException(log, msg);
		}
		
		log.addLine("Working Parameters");
		log.addLine("Output Dir:\t" + args[0]);
		log.addLine("Stats String:\t" + args[1]);
		log.addLine("Window Range:\t" + args[2] + "\n");
	}
	
	private static String[] DEFAULT = {"i", "h", "x", "d", "f"};
	
	private File out_dir;
	private File stats_dir;
	private String[] stat_str;
	
	private List<WindowStats> all_ws;
	
	private Log log;
	
	
	public Combiner(String out_dir_name, Log log) throws IllegalInputException {
		
		this.log = log;
		
		out_dir = new File(out_dir_name);
		if(!out_dir.isDirectory()) {
			String msg = "Error: Output directory path does not exist";
			throw new IllegalInputException(log, msg);
		}
		
		stats_dir = new File(out_dir.getAbsolutePath() + File.separator + "stats_files");
		if(!stats_dir.isDirectory()) {
			String msg = "Error: Stat files directory path does not exist";
			throw new IllegalInputException(log, msg);
		}
		
		all_ws = new LinkedList<WindowStats>();
		stat_str = new String[0];
	}
	
	//combines only some/all of the stats (smarter method)
	public void combineStats(String str, String range) throws Exception {
		
		stat_str = createStatString(str);
		int dwn_rng = getDownRng(range);
		int up_rng = getUpRng(range);
		
		log.addLine("Reading in files");
		importData(dwn_rng, up_rng);
	}
	
	//combines all of the stats
	public void combineStats() {
		
		stat_str = DEFAULT;
	}
	
	//writes the stats to file (can be called from Analyzer too!)
	//NOTE: dont create the output file structure until this method is called
	public void writeStats() {
		
		boolean i,h,x,d,f = false;
		if(ssContains("i"))
			i = true;
		if(ssContains("h"))
			h = true;
		if(ssContains("x"))
			x = true;
		if(ssContains("d"))
			d = true;
		if(ssContains("f"))
			f = true;
	}
	
	public List<WindowStats> getAllStats() {
		return all_ws;
	}
	
	private void importData(int dwn_rng, int up_rng) {
		
		
		
		String[] all_files = stats_dir.list();
		
		for(int i = dwn_rng; i <= up_rng; i++) {
			
			File win_file = getWindowFile(i, all_files);
			if(win_file != null && win_file.exists()) {
				
				int st_pos = getStart(win_file.getName());
				int end_pos = getEnd(win_file.getName());
				
				
				
			}
		}
	}
	
	private int getEnd(String name) {
		
		int st_indx = name.indexOf("-e");
		int end_indx = name.indexOf(".tsv");
		
		System.out.println(name.substring((st_indx + 2), end_indx));
		
		return Integer.parseInt(name.substring((st_indx + 2), end_indx));
	}
	
	private int getStart(String name) {
		
		int st_indx = name.indexOf("_s");
		int end_indx = name.indexOf("-e");
		
		System.out.println(Integer.parseInt(name.substring((st_indx + 2), end_indx)));
		
		return Integer.parseInt(name.substring((st_indx + 2), end_indx));
	}
	
	private File getWindowFile(int win_num, String[] all_files) {
		
		for(int i = 0; i < all_files.length; i++) {
			
			String name = all_files[i];
			if(name.contains("win" + win_num + "_")) 
				return new File(name);
		}
		
		log.addLine("\tWarning: Could not find window number " + win_num);
		return null;
	}
	
	private boolean ssContains(String stat) {
		
		for(int i = 0; i < stat_str.length; i++) {
			if(stat_str[i] == stat)
				return true;
		}
		
		return false;
	}
	
	private int getUpRng(String range) throws IllegalInputException {
		
		String[] rng = range.split("-");
		
		try {
			return Integer.parseInt(rng[1]);
			
		} catch(NumberFormatException e) {
			String msg = "Error: Upper range bound argument " + range + " is invalid";
			throw new IllegalInputException(log, msg);
		}
	}
	
	private int getDownRng(String range) throws IllegalInputException {
		
		String[] rng = range.split("-");
		
		try {
			if(rng.length != 2)
				throw new Exception();
			
			return Integer.parseInt(rng[0]);
			
		} catch(Exception e) {
			String msg = "Error: Lower range bound argument " + range + " is invalid";
			throw new IllegalInputException(log, msg);
		}
	}
	
	private String[] createStatString(String str) throws IllegalInputException {
		
		String[] fin_str = str.split(":");
		
		for(int i = 0; i < fin_str.length; i++) {
			
			String s = fin_str[i];
			if(!s.equals("i") && !s.equals("h") && !s.equals("x") && !s.equals("d") && !s.equals("f")) {
				String msg = "Error: Input stat string has invalid characters or values";
				throw new IllegalInputException(log, msg);
			}
		}
		
		return fin_str;
	}

}
