package analysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedList;
import java.util.List;

import errors.FileParsingException;
import errors.IllegalInputException;
import tools.*;

//WindowCombiner
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
			 c.combineWindows(args[1], args[2]);
			 c.writeStats();
			 
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
	
	private static String[] DEFAULT = {"i", "x", "h", "d", "f"};
	
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
	public void combineWindows(String str, String range) throws Exception {
		
		stat_str = createStatString(str);
		int dwn_rng = getDownRng(range);
		int up_rng = getUpRng(range);
		
		log.addLine("Reading in files");
		importData(dwn_rng, up_rng);
	}
	
	//combines all of the stats
	public void combineWindows() {
		
		stat_str = DEFAULT;
		log.addLine("Reading in files");
		importData();
		
	}
	
	//writes the stats to file (can be called from Analyzer too!)
	//NOTE: dont create the output file structure until this method is called
	public void writeStats() throws FileParsingException {
		
		try {
			out_dir = new File(out_dir.getAbsolutePath() + File.separator + "final_out");
			out_dir.mkdir();
			
			File out_file = new File(out_dir.getAbsoluteFile() + File.separator 
					+ "combined_windows.tsv");
			int num = 1;
			while(out_file.exists()) {
				out_file = new File(out_dir.getAbsoluteFile() + File.separator 
						+ "combined_windows" + num + ".tsv");
				num++;
			}
			out_file.createNewFile();
			
			if(stat_str.equals(DEFAULT))
				simplePrint(out_file);
			else 
				specificPrint(out_file);
			
		} catch(IOException e) {
			String msg = "Error: There was a problem with printing to the output file in " + out_dir.getAbsolutePath();
			throw new FileParsingException(log, msg);
		}
		
	}
	
	public List<WindowStats> getAllStats() {
		return all_ws;
	}
	
	private void specificPrint(File out_file) throws FileNotFoundException {
		
		boolean i,x,h,d,f;
		i = x = h = d = f = false;
		
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
		
		PrintWriter pw = new PrintWriter(out_file);
		pw.print("snp_id\tposition");
		if(i)
			pw.print("\tiHS");
		if(x)
			pw.print("\tXPEHH");
		if(h)
			pw.print("\tiHH");
		if(d)
			pw.print("\tDAF");
		if(f)
			pw.print("\tFst");
		pw.print("\n");
		
		for(int j = 0; j < all_ws.size(); j++) {
			
			WindowStats ws = all_ws.get(j);
			List<SNP> all_snps = ws.getAllSNPs();
			for(int k = 0; k < all_snps.size(); k++) {
				
				SNP snp = all_snps.get(k);
				pw.write(snp.getSnpID() + "\t" + snp.getPosition());
				if(i)
					pw.write("\t" + ws.getIhsScore(snp).toString());
				if(x)
					pw.write("\t" + ws.getXpehhScore(snp).toString());
				if(h)
					pw.write("\t" + ws.getIhhScore(snp).toString());
				if(d)
					pw.write("\t" + ws.getDafScore(snp).toString());
				if(f)
					pw.write("\t" + ws.getFstScore(snp).toString());
				pw.write("\n");	
			}
		}
		
		pw.close();
	}
	
	private void simplePrint(File out_file) throws FileNotFoundException {
		
		PrintWriter pw = new PrintWriter(out_file);
		pw.print("snp_id\tposition\tiHS\tXPEHH\tiHH\tDAF\tFst\n");
		
		for(int i = 0; i < all_ws.size(); i++)
			pw.print(all_ws.get(i));
		
		pw.close();
	}
	
	private void importData() {
		
		String[] all_file_names = stats_dir.list();
		
		for(int i = 0; i < all_file_names.length; i++) {
			
			File win_file = new File(stats_dir.getAbsolutePath() 
					+ File.separator + all_file_names[i]);
			addFile(win_file);
		}
	}
	
	private void importData(int dwn_rng, int up_rng) {
		
		String[] all_files = stats_dir.list();
		
		for(int i = dwn_rng; i <= up_rng; i++) {
			
			File win_file = getWindowFile(i, all_files);
			addFile(win_file);
		}
	}
	
	private void addFile(File win_file) {
		
		if(win_file != null && win_file.exists()) {
			
			int st_pos = getStart(win_file.getName());
			int end_pos = getEnd(win_file.getName());
			
			WindowParser wp = new WindowParser(win_file, st_pos, end_pos);
			WindowStats ws = wp.parseWindow();
			all_ws.add(ws);
		}
	}
	
	private int getEnd(String name) {
		
		int st_indx = name.indexOf("-e");
		int end_indx = name.indexOf(".tsv");
		
		return Integer.parseInt(name.substring((st_indx + 2), end_indx));
	}
	
	private int getStart(String name) {
		
		int st_indx = name.indexOf("_s");
		int end_indx = name.indexOf("-e");
		
		return Integer.parseInt(name.substring((st_indx + 2), end_indx));
	}
	
	private File getWindowFile(int win_num, String[] all_file_names) {
		
		for(int i = 0; i < all_file_names.length; i++) {
			
			String name = all_file_names[i];
			if(name.contains("win" + win_num + "_")) 
				return new File(stats_dir.getAbsolutePath() + File.separator + name);
		}
		
		log.addLine("\tWarning: Could not find window number " + win_num);
		return null;
	}
	
	private boolean ssContains(String stat) {
		
		for(int i = 0; i < stat_str.length; i++) {
			if(stat_str[i].equals(stat))
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