package envi;

import envi.SetupDriver;
import errors.IllegalInputException;
import tools.Log;

public class EnviSetup {

	/** Hyper-Parallelized Composite of Multiple Signals (CMS) Java implementation GWAS and Local study version 1.0
	 * This program sets up the environment needed for stats calculation
	 * @author Hayden Smith
	 * 
	 * @param args0		Phased and Ancestral Data Files Directory Path (.phased/.legend -OR- .vcf)
	 * @param args1		Genetic Map File Path
	 * @param args2		Environment Output Directory Path
	 * @param args3		Target Population (CEU, YRI, JPT, CHB)
	 * @param args4		Cross Population (CEU, YRI, JPT, CHB)
	 * @param args5		Outgroup Population (CEU, YRI, JPT, CHB)
	 * @param args6		Chromosome Declaration (i.e. 1-22)
	 * @param args7		Window Size in Mb
	 */
	public static void main(String[] args) {
		
		Log log = new Log(Log.type.envi);
		
		try {
			args = setupArgs(args, log);
		
			SetupDriver dv = new SetupDriver(args, log);
			dv.runSetup();
			
			System.out.println("Setup Finished");
			log.addLine("Your environment has been successfully set up. Go to api "
					+ "for information about how to run stats analysis");
			
		} catch (Exception e) {
			
			System.out.println("CMS Died Prematurely." 
					+ " Check log output for troubleshooting.");
			
			log.addLine("\n\nCMS Died Prematurely. Error in computation.");
			
			e.printStackTrace();
		}
	}
	
	private static String[] setupArgs(String[] args, Log log) 
			throws IllegalInputException{
	
		if(args.length != 8) {
			String msg = "Error: Parameter length incorrect";
			throw new IllegalInputException(log, msg);
		}
	
		log.addLine("Working Parameters");
		log.addLine("Data Dir:\t\t" + args[0]);
		log.addLine("Map File:\t\t" + args[1]);
		log.addLine("Envi Output Dir:\t" + args[2]);
		log.addLine("Target Pop:\t\t" + args[3]);
		log.addLine("Cross Pop:\t\t" + args[4]);
		log.addLine("Outgroup Pop:\t\t" + args[5]);
		log.addLine("Chr Range:\t\t" + args[6]);
		log.addLine("Window Size:\t\t" + args[7] + "Mb");
		
		return args;
	}
}
