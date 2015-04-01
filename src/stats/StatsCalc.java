package stats;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.List;

import calc.*;
import tools.*;

public class StatsCalc {

	/** Hyper-Parallelized Composite of Multiple Signals (CMS) Java implementation GWAS and Local study version 1.0
	 * This program calculated the stats that are used in CMS analysis
	 * @author Hayden Smith
	 * 
	 * @param args0		out dir
	 * @param args1		chr
	 * @param args2		window number
	 */
	public static void main(String[] args) {

		Log log = new Log(Log.type.stat);
		
		File out_dir = new File(args[0]);
		if(!out_dir.exists() && !out_dir.isDirectory()) {
			//TODO: throw some error
		}
		
		int chr = Integer.parseInt(args[1]);
		int win_num = Integer.parseInt(args[2]);
		
		StatsCalc sc = new StatsCalc(out_dir, chr, win_num, log);
		sc.runStats();
		
		log.addLine("Run Complete!");
	}
	
	private static int WAIT_TIME = 50;
	
	private int chr;
	private int win_num;
	
	private File out_dir;
	private File envi_dir;
	private File win_dir;
	
	private iHS i;
	private iHH h;
	private XPEHH x;
	private DAF d;
	private Fst f;
	
	private Log log;
	
	public StatsCalc(File out_dir, int chr, int win_num, Log log) {
		
		this.chr = chr;
		this.win_num = win_num;
		this.out_dir = out_dir;
		this.log = log;
	}
	
	public void runStats() {
		
		setupFiles();
		createCalculators();
		doCalculations();
		
		//run the tests
		//write the stats output file
		writeOutput();
	}
	
	private void writeOutput() {
		
		//TODO:
		i.getStats();
		i.getSNPs();
		h.getStats();
		h.getSNPs();
		x.getStats();
		x.getSNPs();
		d.getStats();
		d.getSNPs();
		f.getStats();
		f.getSNPs();
	}
	
	private void doCalculations() {
		
		Object lock = new Object();
		
		try {
		
			StatsThread i_thrd = new StatsThread(i, lock);
			Thread.sleep(WAIT_TIME);
			i_thrd.start();
			
			StatsThread h_thrd = new StatsThread(h, lock);
			Thread.sleep(WAIT_TIME);
			h_thrd.start();
			
			StatsThread x_thrd = new StatsThread(x, lock);
			Thread.sleep(WAIT_TIME);
			x_thrd.start();
			
			StatsThread d_thrd = new StatsThread(d, lock);
			Thread.sleep(WAIT_TIME);
			d_thrd.start();
			
			StatsThread f_thrd = new StatsThread(f, lock); 
			Thread.sleep(WAIT_TIME);
			f_thrd.start();
			
			synchronize(i_thrd, h_thrd, x_thrd, d_thrd, f_thrd);
			
			i = (iHS) i_thrd.getTest();
			h = (iHH) h_thrd.getTest();
			x = (XPEHH) x_thrd.getTest();
			d = (DAF) d_thrd.getTest();
			f = (Fst) f_thrd.getTest();
		
		} catch (InterruptedException e) { //TODO: throw new threading error
			e.printStackTrace();
		}
		
		i.logRunStats();
		h.logRunStats();
		x.logRunStats();
		d.logRunStats();
		f.logRunStats();
	}
	
	private void synchronize(StatsThread i_thrd, 
								StatsThread h_thrd, 
								StatsThread x_thrd, 
								StatsThread d_thrd, 
								StatsThread f_thrd) {

		for(;;) {
			if(i_thrd.isFinished()
					&& h_thrd.isFinished()
					&& x_thrd.isFinished()
					&& d_thrd.isFinished()
					&& f_thrd.isFinished())
				break;
			else
				continue;
		}
	}
	
	@SuppressWarnings("unchecked")
	private void createCalculators() {
		
		try { 
			String path = "";
			
			path = getEnviFileName("target_pop_wins.bin");
			List<Window> tp_wins = (List<Window>) getObject(path);
			
			path = getEnviFileName("target_pop_indv.bin");
			Individual[] tp_indv = (Individual[]) getObject(path);
			
			path = getEnviFileName("targetXcross_wins.bin");
			List<Window> txin_wins = (List<Window>) getObject(path);
			
			path = getEnviFileName("targetXcross_indv.bin");
			Individual[] tp_inx_indv = (Individual[]) getObject(path);
			
			path = getEnviFileName("crossXtarget_indv.bin");
			Individual[] xp_int_indv = (Individual[]) getObject(path);
			
			path = getEnviFileName("crossXout_wins.bin");
			List<Window> xoin_wins = (List<Window>) getObject(path);
			
			path = getEnviFileName("crossXout_indv.bin");
			Individual[] xp_ino_indv = (Individual[]) getObject(path);
			
			path = getEnviFileName("outXcross_indv.bin");
			Individual[] op_inx_indv = (Individual[]) getObject(path);
			
			path = getEnviFileName("anc_types.bin");
			List<Window> anc_types = (List<Window>) getObject(path);
			
			path = getEnviFileName("genetic_map.bin");
			GeneticMap gm = (GeneticMap) getObject(path);
			
			path = getTargetWindowFileName();
			Window tp_win = (Window) getObject(path);
			
			path = getTargetXCrossWindowFileName();
			Window txin_win = (Window) getObject(path);
			
			i = new iHS(log, tp_win, tp_indv, anc_types, tp_wins, gm);
			h = new iHH(log, tp_win, tp_indv, anc_types, tp_wins, gm);
			x = new XPEHH(log, txin_win, txin_wins, tp_inx_indv, xp_int_indv, gm);
			d = new DAF(log, tp_win, tp_indv, xoin_wins, xp_ino_indv, op_inx_indv, anc_types);
			f = new Fst(log, txin_win, tp_inx_indv, xp_int_indv, op_inx_indv);
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClassCastException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	@SuppressWarnings("resource")
	private Object getObject(String path) throws FileNotFoundException, IOException, ClassNotFoundException {
		
		File file = new File(path);
		if(!file.exists()) {
			//TODO: throw some error
		}
		
		ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(file)));
		return ois.readObject();
	}
	
	private String getTargetXCrossWindowFileName() throws FileNotFoundException {
		
		String[] all_files = win_dir.list();
		for(int i = 0; i < all_files.length; i++) {
			
			String file_name = all_files[i];
			if(file_name.contains("chr" + chr)
					&& file_name.contains("win" + win_num)
					&& file_name.contains("x"))
				return file_name;
		}
		
		//TODO: change to make better??? be specific about what file
		throw new FileNotFoundException();
	}
	
	private String getTargetWindowFileName() throws FileNotFoundException {
		
		String[] all_files = win_dir.list();
		for(int i = 0; i < all_files.length; i++) {
			
			String file_name = all_files[i];
			if(file_name.contains("chr" + chr)
					&& file_name.contains("win" + win_num)
					&& !file_name.contains("x"))
				return file_name;
		}
		
		//TODO: change to make better??? be specific about what file
		throw new FileNotFoundException();
	}
	
	private String getEnviFileName(String name) {
		return envi_dir.getAbsolutePath() + File.separator + name;
	}
	
	private void setupFiles() {
		
		envi_dir = new File(out_dir.getAbsoluteFile() + File.separator + "envi_files" + File.separator + "envi_var");
		if(!envi_dir.exists() && !envi_dir.isDirectory()) {
			//TODO: throw some error
		}
		
		win_dir = new File(out_dir.getAbsoluteFile() + File.separator + "envi_files" + File.separator + "all_wins");
		if(!win_dir.exists() && !win_dir.isDirectory()) {
			//TODO: throw some error
		}
		
		out_dir = new File(out_dir.getAbsolutePath() + File.separator + "stats_files");
		if(!out_dir.exists())
			out_dir.mkdirs();
		
		
	}
	
	

}

class StatsThread extends Thread {

	private final Object lock;
	
	private Thread thrd;
	private HaplotypeTests tst;
	
	volatile private boolean finished;//volatile says that some other thread could change this value
	
	StatsThread(HaplotypeTests tst, Object lock) {
		this.tst = tst;
		this.lock = lock;
		
		finished = false;
		
		thrd = new Thread(this);
//		thrd.start();
	}
	
	public void start() {
		thrd.start();
	}
	
	@Override
	public void run() {
		
		tst.runStat();	
		
		synchronized(lock) {
			finished = true;
		}
		
		thrd.interrupt();
	}
	
	public HaplotypeTests getTest() {
		return tst;
	}
	
	public boolean isFinished() {
		
		synchronized(lock) {
			return finished;
		}
	}
	
}
