package analysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.special.Erf;

import errors.IllegalInputException;
import tools.Log;
import tools.SNP;
import tools.WindowStats;

public class SignificanceAnalyzer {

	/**
	 * Gets all the data from all the stats output files
	 * If directly called from main it will print the data (any combination of tests) into one big file
	 * If called from Analyzer it will return List<WindowStats> for further analysis
	 * 
	 * FOR WINDOW COMBINING AND NO STATS ANALYSIS
	 * @param args0		window combining flag ("combiner")
	 * @param args1		out dir
	 * @param args2		chr
	 * @param args3		stats string ("i:h:x:d:f" or any combination of that) (optional)
	 * 
	 * FOR STATS ANALYSIS ONLY
	 * @param args0		out dir
	 * @param args1		chr
	 * @param args2		p-value (optional)
	 */
	public static void main(String[] args) {
		
		Log log = new Log(Log.type.analysis);
		System.out.println("Starting the final phase of SelecT");
		
		String type = args[0];
		
		try {
			
			int chr = -1;
			try {
				
				if(type.equals("combine"))
					chr = Integer.parseInt(args[2]);
				else
					chr = Integer.parseInt(args[1]);
				
				if(chr < 1 || chr > 22) {
					String msg = "Error: Chromosome declaration is out of bounds";
					throw new IllegalInputException(log, msg);
				}
				
			} catch(NumberFormatException e) {
				String msg = "Error: Invalid chromosome format";
				throw new IllegalInputException(log, msg);
			}
			
			if(type.equals("combine")) {
				
				checkCombineArgs(args, log);
				File wrk_dir = new File(args[1]);
				if(!wrk_dir.isDirectory()) {
					String msg = "Error: Output directory path does not exist";
					throw new IllegalInputException(log, msg);
				}
				
				Combiner c = new Combiner(wrk_dir, chr, log);
				
				if(args.length == 3)
					c.combineWindows();
				else
					c.combineWindows(args[3]);
				
				c.writeStats();
			}
			else {
				
				checkStatsAnalArgs(args, log);
				File wrk_dir = new File(args[0]);
				if(!wrk_dir.isDirectory()) {
					String msg = "Error: Output directory path does not exist";
					throw new IllegalInputException(log, msg);
				}
				
				double p_val = 0.01;//TODO: allow for p-value exploration (default = 0.01)
				boolean write_all_win_to_one_file = false;
				boolean run_normalization = false;
				boolean ignore_mop = false;
				boolean ignore_pop = false;
				try {
					p_val = Double.parseDouble(args[2]);
					
					if(args[3].equals("T"))
						write_all_win_to_one_file = true;
					else if(args[3].equals("F"))
						write_all_win_to_one_file = false;
					else {
						String msg = "The fourth argument is in illegal format";
						throw new IllegalInputException(log, msg);
					}
					
					if(args[4].equals("T"))
						run_normalization = true;
					else if(args[4].equals("F"))
						run_normalization = false;
					else {
						String msg = "The fifth argument is in illegal format";
						throw new IllegalInputException(log, msg);
					}
					
					if(args[5].equals("T"))
						ignore_mop = true;
					else if(args[5].equals("F"))
						ignore_mop = false;
					else {
						String msg = "The sixth argument is in illegal format";
						throw new IllegalInputException(log, msg);
					}
					
					if(args[6].equals("T"))
						ignore_pop = true;
					else if(args[6].equals("F"))
						ignore_pop = false;
					else {
						String msg = "The seventh argument is in illegal format";
						throw new IllegalInputException(log, msg);
					}
					
					if(ignore_mop && ignore_pop) {//TODO: make this case an either or situation
						String msg = "Error: Ignoring both MoP and PoP scores is an invalid request for analysis";
						throw new IllegalInputException(log, msg);
					}
					
				} catch(NumberFormatException e) {
					String msg = "Error: Invalid p-value input";
					throw new IllegalInputException(log, msg);
				}
				
				Combiner c = new Combiner(wrk_dir, chr, log);
				c.combineAnalysisData();
				
				if(write_all_win_to_one_file)
					c.writeStats();
				
				System.out.println("\nStarting Analysis");
				SignificanceAnalyzer sa = new SignificanceAnalyzer(log, p_val, c.getAllStats(), wrk_dir);
				sa.findSignificantSNPs(run_normalization, ignore_mop, ignore_pop);
			}
			 
		} catch (Exception e) {
			
			System.out.println("CMS Died Prematurely." 
					+ " Check log output for troubleshooting.");
			
			log.addLine("\n\nCMS Died Prematurely. Error in computation.");
			
			e.printStackTrace();
		}	
		
		System.out.println("\n\nSelecT significance analysis finished!");
	}
	
	private static void checkStatsAnalArgs(String[] args, Log log) 
			throws IllegalInputException {
		
		if(args.length != 7) {
			String msg = "Error: Parameter length incorrect";
			throw new IllegalInputException(log, msg);
		}
		
		log.addLine("Working Parameters");
		log.addLine("Output Dir:\t" + args[0]);
		log.addLine("Chromosome:\t" + args[1]);
		log.addLine("p-Value:\t" + args[2]);
	}
	
	private static void checkCombineArgs(String[] args, Log log) 
			throws IllegalInputException {
		
		if(args.length < 3 && args.length > 4) {
			String msg = "Error: Parameter length incorrect";
			throw new IllegalInputException(log, msg);
		}
		
		log.addLine("Working Parameters");
		log.addLine("Output Dir:\t" + args[1]);
		log.addLine("Chromosome:\t" + args[2]);
		if(args.length == 4)
			log.addLine("Stats String:\t" + args[3] + "\n");
	}
	
	
	private double sig_score;
	private File out_file;
	private List<WindowStats> all_ws;
	
	private Log log;
	
	public SignificanceAnalyzer(Log log, double sig_pval, List<WindowStats> all_ws, File wrk_dir) {
		
		this.log = log;
		this.all_ws = all_ws;
		
		sig_score = pToZ(sig_pval);
		
		System.out.println("p-value:\t\t" + sig_pval);
		System.out.println("Significant Score:\t" + sig_score);
		
		wrk_dir = new File(wrk_dir.getAbsolutePath() + File.separator + "final_out");
		if(!wrk_dir.exists())
			wrk_dir.mkdir();
		out_file = new File(wrk_dir.getAbsoluteFile() + File.separator + "significant_loci.tsv");
		int num = 1;
		while(out_file.exists()) {
			out_file = new File(wrk_dir.getAbsoluteFile() + File.separator 
					+ "significant_loci" + num + ".tsv");
			num++;
		}
	}
	
	public void findSignificantSNPs(boolean run_normalization, boolean ignore_mop, boolean ignore_pop) throws IllegalInputException {
		
		if(run_normalization) 
			all_ws = normalizeAllWindows(all_ws);
		
		try {
			
			System.out.println("Extracting significant loci");
			
			PrintWriter pw = new PrintWriter(out_file);
			pw.print("snp_id\tposition\tiHS\tXPEHH\tiHH\tdDAF\tDAF\tFst"//\tTajD\tNew
					+ "\tunstd_PoP\tunstd_MoP\twin_PoP\twin_MoP\n");
			
			for(int i = 0; i < all_ws.size(); i++) {
				
				WindowStats ws = all_ws.get(i);
				List<SNP> ws_snps = ws.getAllSNPs();
				
				for(int j = 0; j < ws_snps.size(); j++) {
					
					if(ignore_mop) {
						if(ws.getStdPopScore(ws_snps.get(j)) >= sig_score)
							pw.print(ws.printSNP(ws_snps.get(j)));
					} else if(ignore_pop) {
						if(ws.getStdMopScore(ws_snps.get(j)) >= sig_score)
							pw.print(ws.printSNP(ws_snps.get(j)));
					} else {
						if(ws.getStdPopScore(ws_snps.get(j)) >= sig_score
								&& ws.getStdMopScore(ws_snps.get(j)) >= sig_score)
							pw.print(ws.printSNP(ws_snps.get(j)));
					}
				}
			}
			
			pw.close();
			
		} catch (FileNotFoundException e) {
			String msg = "Error: Could not create output file for final analysis";
			throw new IllegalInputException(log, msg);
		}
	}
	
	private List<WindowStats> normalizeAllWindows(List<WindowStats> all_stats) {
		
		System.out.println("Running Normalization");
		
		WindowStats comb_ws = new WindowStats(all_stats.get(0).getStPos(), 
				all_stats.get(all_stats.size()-1).getEndPos());
		
		for(int i = 0; i < all_stats.size(); i++) {
			
			WindowStats ws = all_stats.get(i);
			comb_ws.addIHS(ws.getIHSstats(), ws.getIHSsnps());
			comb_ws.addXPEHH(ws.getXPEHHstats(), ws.getXPEHHsnps());
			comb_ws.addIHH(ws.getIHHstats(), ws.getIHHsnps());
			comb_ws.addDDAF(ws.getDDAFstats(), ws.getDDAFsnps());
			comb_ws.addDAF(ws.getDAFstats(), ws.getDAFsnps());
			comb_ws.addFST(ws.getFSTstats(), ws.getFSTsnps());
			//comb_ws.addTAJD(ws.getTAJDstats(), ws.getTAJDsnps());
			//comb_ws.addNEW(ws.getNEWstats(), ws.getNEWsnps());
			
			comb_ws.addUnstdPoP(ws.getUnstdPoP());
			comb_ws.addUnstdMoP(ws.getUnstdMoP());
		}
		
		comb_ws.normalizeUnstdCompositeScores();
		List<WindowStats> comb_ws_list = new ArrayList<WindowStats>();
		comb_ws_list.add(comb_ws);
		
		return comb_ws_list;
	}
	
	private double pToZ(double p) {
	    double z = Math.sqrt(2) * Erf.erfcInv(2*p);
	    return(z);
	}
}
