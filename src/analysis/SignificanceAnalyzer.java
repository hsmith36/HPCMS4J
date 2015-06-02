package analysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;

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
	 * @param args0		SelecT workspace
	 * @param args1		chr
	 * @param args2		p-value (optional)
	 */
	public static void main(String[] args) {
		
		System.out.println("Starting the final phase of SelecT");
		
		Log log = new Log(Log.type.analysis);
		HashMap<String, Object> arg_map = setupArgs(args);
		
//		String type = args[0];
		
		try {
			
//			int chr = -1;
//			try {
//				
//				if(type.equals("combine"))
//					chr = Integer.parseInt(args[2]);
//				else
//					chr = Integer.parseInt(args[1]);
//				
//				if(chr < 1 || chr > 22) {
//					String msg = "Error: Chromosome declaration is out of bounds";
//					throw new IllegalInputException(log, msg);
//				}
//				
//			} catch(NumberFormatException e) {
//				String msg = "Error: Invalid chromosome format";
//				throw new IllegalInputException(log, msg);
//			}
//			
//			if(type.equals("combine")) {
//				
//				checkCombineArgs(args, log);
//				File wrk_dir = new File(args[1]);
//				if(!wrk_dir.isDirectory()) {
//					String msg = "Error: Output directory path does not exist";
//					throw new IllegalInputException(log, msg);
//				}
//				
//				Combiner c = new Combiner(wrk_dir, chr, log);
//				
//				if(args.length == 3)
//					c.combineWindows();
//				else
//					c.combineWindows(args[3]);
//				
//				c.writeStats();
//			}
//			else {
//				
//				checkStatsAnalArgs(args, log);
//				File wrk_dir = new File(args[0]);
//				if(!wrk_dir.isDirectory()) {
//					String msg = "Error: Output directory path does not exist";
//					throw new IllegalInputException(log, msg);
//				}
//				
//				double p_val = 0.01;
//				boolean write_all_win_to_one_file = false;
//				boolean run_normalization = false;
//				boolean ignore_mop = false;
//				boolean ignore_pop = false;
//				try {
//					p_val = Double.parseDouble(args[2]);
//					
//					if(args[3].equals("T"))
//						write_all_win_to_one_file = true;
//					else if(args[3].equals("F"))
//						write_all_win_to_one_file = false;
//					else {
//						String msg = "The fourth argument is in illegal format";
//						throw new IllegalInputException(log, msg);
//					}
//					
//					if(args[4].equals("T"))
//						run_normalization = true;
//					else if(args[4].equals("F"))
//						run_normalization = false;
//					else {
//						String msg = "The fifth argument is in illegal format";
//						throw new IllegalInputException(log, msg);
//					}
//					
//					if(args[5].equals("T"))
//						ignore_mop = true;
//					else if(args[5].equals("F"))
//						ignore_mop = false;
//					else {
//						String msg = "The sixth argument is in illegal format";
//						throw new IllegalInputException(log, msg);
//					}
//					
//					if(args[6].equals("T"))
//						ignore_pop = true;
//					else if(args[6].equals("F"))
//						ignore_pop = false;
//					else {
//						String msg = "The seventh argument is in illegal format";
//						throw new IllegalInputException(log, msg);
//					}
//					
////					if(ignore_mop && ignore_pop) {
////						String msg = "Error: Ignoring both MoP and PoP scores is an invalid request for analysis";
////						throw new IllegalInputException(log, msg);
////					}
//					
//				} catch(NumberFormatException e) {
//					String msg = "Error: Invalid p-value input";
//					throw new IllegalInputException(log, msg);
//				}
//				
//				Combiner c = new Combiner(wrk_dir, chr, log);
//				c.combineAnalysisData();
//				
//				if(write_all_win_to_one_file)
//					c.writeStats();
//				
//				System.out.println("\nStarting Analysis");
//				SignificanceAnalyzer sa = new SignificanceAnalyzer(log, p_val, c.getAllStats(), wrk_dir);
//				sa.findSignificantSNPs(run_normalization, ignore_mop, ignore_pop);
//			}
			
			//if arg_map.getType == combiner 
			//	create Combiner
			//	run Combiner on optional string (if not present run combineWindow())
			//else
			//	create Combiner
			//	run generic combineWindows()
			//	create SignificanceAnalyzer
			//	run analyzer
			
			if((Boolean) arg_map.get("combine_only")) {
				
				Combiner c = new Combiner(arg_map, log);
				if(arg_map.get("combine_filter").equals(".:default:."))
					c.combineWindows();
				else
					c.combineWindows((String) arg_map.get("combine_filter"));
				
				c.writeStats();
			}
			else {
				
				Combiner c = new Combiner(arg_map, log);
				c.combineAnalysisData();
				
				if((Boolean) arg_map.get("write_combine"))
					c.writeStats();
				
				System.out.println("\nStarting Analysis");
				boolean run_normalization = (Boolean) arg_map.get("run_norm");
				boolean ignore_mop = (Boolean) arg_map.get("ignore_mop");
				boolean ignore_pop = (Boolean) arg_map.get("ignore_pop");
				
				SignificanceAnalyzer sa = new SignificanceAnalyzer(arg_map, c.getAllStats(), log);
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
	
	private static HashMap<String, Object> setupArgs(String[] args) {
		
		ArgumentParser parser = ArgumentParsers.newArgumentParser("SignificanceAnalyzer")
				.defaultHelp(true)
                .description("Run Analysis on evolutionary statistics");
		
		//Creating required arguments
		parser.addArgument("wrk_dir").type(Arguments.fileType().verifyIsDirectory()
                .verifyCanRead()).help("SelecT workspace directory");
		
		parser.addArgument("chr").type(Integer.class).choices(Arguments.range(1, 22))
				.help("Chromosome number");
		
		//Creating optional arguments
		parser.addArgument("-co", "--combine_only").action(Arguments.storeTrue())
				.help("When present this flag only runs fist half of analysis where "
						+ "windows are combined into one file");
		
		parser.addArgument("-cf", "--combine_filter").type(String.class).setDefault(".:default:.")
				.help("Uses a specific filter for specific combination of stats in output file "
						+ "and can only be used in conjunction with the --combine_only flag. "
						+ "See api for definition of filter string");
		
		parser.addArgument("-p", "--p_value").type(Double.class).setDefault(0.01)
				.choices(Arguments.range(0.000000000000000000001, 0.99999999999999999))
				.help("Sets the p-value cutoff for significance check on composite scores. "
						+ "If not included, defaults to 0.01");
		
		parser.addArgument("-wc", "--write_combine").action(Arguments.storeTrue())
				.help("During full analysis this flag creates a separate file "
						+ "complete with all unfiltered data");
		
		parser.addArgument("-rn", "--run_norm").action(Arguments.storeTrue())
				.help("Normalizes MoP and PoP probabilities across whole dataset instead of "
						+ "one window before finding significance");
		
		parser.addArgument("-im", "--ignore_mop").action(Arguments.storeTrue())
				.help("Ignores all MoP scores and finds significance based upon PoP only; "
						+ "when both --ignore flags are present significance is determined "
						+ "by finding either MoP or PoP significance");
		
		parser.addArgument("-ip", "--ignore_pop").action(Arguments.storeTrue())
				.help("Ignores all PoP scores and finds significance based upon MoP only; "
						+ "when both --ignore flags are present significance is determined "
						+ "by finding either MoP or PoP significance");
		
		//Parsing user-inputed arguments
		HashMap<String, Object> parsedArgs = new HashMap<String, Object>();
		
		//Checking to make sure input is correct
		try {parser.parseArgs(args, parsedArgs);}
    	catch (ArgumentParserException e) {
    		System.out.println("Fatal error in argument parsing: see log");
            e.printStackTrace();
            
            Log err_log = new Log(Log.type.stat);
            err_log.addLine("Error: Failed to parse arguments"); 
            err_log.addLine("\t*" + e.getMessage());
            err_log.addLine("\t*Go to api for more information or run with -h as first parameter for help");
			err_log.addLine("\t*You will need to redo this entire step--all new data is invalid");
			
			System.exit(0);
        }
		
		return parsedArgs;
	}
	
//	private static void checkStatsAnalArgs(String[] args, Log log) 
//			throws IllegalInputException {
//		
//		if(args.length != 7) {
//			String msg = "Error: Parameter length incorrect";
//			throw new IllegalInputException(log, msg);
//		}
//		
//		log.addLine("Working Parameters");
//		log.addLine("Output Dir:\t" + args[0]);
//		log.addLine("Chromosome:\t" + args[1]);
//		log.addLine("p-Value:\t" + args[2]);
//	}
//	
//	private static void checkCombineArgs(String[] args, Log log) 
//			throws IllegalInputException {
//		
//		if(args.length < 3 && args.length > 4) {
//			String msg = "Error: Parameter length incorrect";
//			throw new IllegalInputException(log, msg);
//		}
//		
//		log.addLine("Working Parameters");
//		log.addLine("Output Dir:\t" + args[1]);
//		log.addLine("Chromosome:\t" + args[2]);
//		if(args.length == 4)
//			log.addLine("Stats String:\t" + args[3] + "\n");
//	}
	
	
	private double sig_score;
	private File out_file;
	private List<WindowStats> all_ws;
	
	private Log log;
	
//	public SignificanceAnalyzer(Log log, double sig_pval, List<WindowStats> all_ws, File wrk_dir) {
//		
//		this.log = log;
//		this.all_ws = all_ws;
//		
//		sig_score = pToZ(sig_pval);
//		
//		System.out.println("p-value:\t\t" + sig_pval);
//		System.out.println("Significant Score:\t" + sig_score);
//		
//		wrk_dir = new File(wrk_dir.getAbsolutePath() + File.separator + "final_out");
//		if(!wrk_dir.exists())
//			wrk_dir.mkdir();
//		out_file = new File(wrk_dir.getAbsoluteFile() + File.separator + "significant_loci.tsv");
//		int num = 1;
//		while(out_file.exists()) {
//			out_file = new File(wrk_dir.getAbsoluteFile() + File.separator 
//					+ "significant_loci" + num + ".tsv");
//			num++;
//		}
//	}
	
	public SignificanceAnalyzer(HashMap<String, Object> arg_map, List<WindowStats> all_ws, Log log) {
		
		this.log = log;
		this.all_ws = all_ws;
		
		sig_score = pToZ((Double) arg_map.get("p_value"));
		System.out.println("Using p-value of " + (Double) arg_map.get("p_value") 
				+ " and significance score of " + sig_score);
		
		File wrk_dir = (File) arg_map.get("wrk_dir");
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
			
			out_file.createNewFile();
			PrintWriter pw = new PrintWriter(out_file);
			pw.print("snp_id\tposition\tiHS\tXPEHH\tiHH\tdDAF\tDAF\tFst"//\tTajD\tNew
					+ "\tunstd_PoP\tunstd_MoP\twin_PoP\twin_MoP\n");
			pw.flush();
			
			System.out.println("here");
			
			for(int i = 0; i < all_ws.size(); i++) {
				
				WindowStats ws = all_ws.get(i);
				List<SNP> ws_snps = ws.getAllSNPs();
				
				for(int j = 0; j < ws_snps.size(); j++) {
					
					SNP s = ws_snps.get(j);
					
					if(ignore_mop && ignore_pop) {
						if(ws.getStdPopScore(s) >= sig_score
								|| ws.getStdMopScore(s) >= sig_score) {
							pw.print(ws.printSNP(s));
							pw.flush();
						}
					} else if(ignore_mop) {
						if(ws.getStdPopScore(s) >= sig_score) {
							pw.print(ws.printSNP(s));
							pw.flush();
						}
					} else if(ignore_pop) {
						if(ws.getStdMopScore(s) >= sig_score) {
							pw.print(ws.printSNP(s));
							pw.flush();
						}
					} else {
						if(ws.getStdPopScore(s) >= sig_score
								&& ws.getStdMopScore(s) >= sig_score) {
							pw.print(ws.printSNP(s));
							pw.flush();
						}
					}
				}
			}
			
			pw.close();
			
		} catch (FileNotFoundException e) {
			String msg = "Error: Could not create output file for final analysis";
			throw new IllegalInputException(log, msg);
		} catch (IOException e) {
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
//			comb_ws.addIHS(ws.getIHSstats(), ws.getIHSsnps());
//			comb_ws.addXPEHH(ws.getXPEHHstats(), ws.getXPEHHsnps());
//			comb_ws.addIHH(ws.getIHHstats(), ws.getIHHsnps());
//			comb_ws.addDDAF(ws.getDDAFstats(), ws.getDDAFsnps());
//			comb_ws.addDAF(ws.getDAFstats(), ws.getDAFsnps());
//			comb_ws.addFST(ws.getFSTstats(), ws.getFSTsnps());
			//comb_ws.addTAJD(ws.getTAJDstats(), ws.getTAJDsnps());
			//comb_ws.addNEW(ws.getNEWstats(), ws.getNEWsnps());
			
			comb_ws.addUnstdPoP(ws.getUnstdPoP());
			comb_ws.addUnstdMoP(ws.getUnstdMoP());
		}
		
		comb_ws.normalizeUnstdCompositeScores();
//		List<WindowStats> comb_ws_list = new ArrayList<WindowStats>();
		
		TreeMap<SNP, Double> std_pop = comb_ws.getStdPoP();
		Set<SNP> all_pop_snps = std_pop.keySet();
		for(SNP s : all_pop_snps) {
			for(int i = 0; i < all_stats.size(); i++) {
				if(all_stats.get(i).containsSNP(s))
					all_stats.get(i).addStdPopScore(s, std_pop.get(s));
			}
		}
		
		TreeMap<SNP, Double> std_mop = comb_ws.getStdMoP();
		Set<SNP> all_mop_snps = std_pop.keySet();
		for(SNP s : all_mop_snps) {
			for(int i = 0; i < all_stats.size(); i++) {
				if(all_stats.get(i).containsSNP(s))
					all_stats.get(i).addStdMopScore(s, std_mop.get(s));
			}
		}
		
		
//		comb_ws_list.add(comb_ws);
		
		return all_stats;
	}
	
	private double pToZ(double p) {
	    double z = Math.sqrt(2) * Erf.erfcInv(2*p);
	    return(z);
	}
}
