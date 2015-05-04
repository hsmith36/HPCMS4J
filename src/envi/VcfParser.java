package envi;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import errors.FileParsingException;
import tools.Individual;
import tools.Window;
import log.Log;



public class VcfParser {
	
	private String file_path;
	private int chr; 
	private List<Window> windows;
	private List<Window> ancestral;
	private Individual[] individuals;
	private Log log;
	
	public VcfParser()
	{
		file_path = null; 
		chr = 0; 
		windows = null; 
		ancestral = null; 
		individuals = null;
		log = null; 
	}
	
	public VCFParser(String file_path, int chr, Log log)
	{
		this.file_path = file_path; 
		this.chr = chr; 
		this.log = log;
		
		individuals = new Individual[0]; // Place holder
		windows = new ArrayList<Window>();
		ancestral = new ArrayList<Window>();
	}
	
	public void parseVCF(int wind_len, boolean ancestral_parser) throws FileParsingException
	{
		try (BufferedReader br = new BufferedReader(new FileReader(file_path))) 
		{
		    String line;
		    while ((line = br.readLine()) != null) //get to header line
		    {
		       if (line.charAt(1) != '#') // line with header titles
		       {
		    	   List<String> indiv_ids = new ArrayList<String>(Arrays.asList(line.split("\t")));
		    	   individuals = new Individual[indiv_ids.size()-9];
		    	   for (int i = 9; i < indiv_ids.size(); i++)
		    	   {
		    		   individuals[i-9] = new Individual(Integer.parseInt(indiv_ids.get(i).substring(2)), chr); // You'll have to change chr type to int in the Individual class for this to work
		    	   }
		    	   break;
		       }   
		    }
		    
		    Window window = new Window();
		    Window anc_window = new Window();
		    int start_pos = 0;
		    int end_pos = wind_len;
		    int index = 0;
		    int pos = 0;
		    
		    while ((line = br.readLine()) != null) // parse the rest of the file
		    {
		    	index++; 
		    	List<String> line_list = new ArrayList<String>(Arrays.asList(line.split("\t")));
		    	
		    	// UPDATE WINDOWS
		    	pos = Integer.parseInt(line_list.get(1));
		    	if (pos > (start_pos + wind_len)) // Check to see if it's time to start a new window 
		    	{
		    		if(!window.equals(new Window()))
		    			windows.add(window);
		    		
		    		if (ancestral_parser == true)
		    		{
		    			if(anc_window.getSNPs() != null && anc_window.getSnpListSize() > 0)
			    			ancestral.add(anc_window);
		    		}
		    		
		    		while (pos > (start_pos + wind_len))
		    		{
		    			start_pos += wind_len;
		    			end_pos += wind_len;
		    		}
		    		
		    		window.setEndIndex(index - 1);
		    		window = new Window(start_pos, end_pos, index);
		    		if (ancestral_parser == true)
		    		{
		    			anc_window.setEndIndex(index - 1);
			    		anc_window = new Window(start_pos, end_pos, index);
		    		}
		    	}
		    
		    	//get SNP information & add to current window
		    	window.addSNP(pos, line_list.get(3), line_list.get(4), line_list.get(2));
		    	if (ancestral_parser && validAncestralData(line_list.get(7)))
		    		anc_window.addSNP(pos, ancestralAllele(line_list), "-", line_list.get(2)); //A1 value for ancestral is -
		       
		    	//UPDATE INVIDUALS
		    	for (int i = 9; i < line_list.size(); i++) // 9 is the column number where individual info starts 
		    	{
		    		List<String> alleles = new ArrayList<String>(Arrays.asList(line_list.get(i).split("\\|")));
		    		individuals[i-9].addAlleleToStrand1(alleles.get(0)); 
		    		individuals[i-9].addAlleleToStrand2(alleles.get(1)); 
		    	}
		    }
		    
		    window.setEndIndex(index - 1);
			windows.add(window);
			
			
			
			
		    // PRINT TESTS
		    /*System.out.println("INDIVIDUALS");
			for (int i = 0; i < individuals.length; i++)
			{
				System.out.println(i);
				System.out.println(individuals[i].toString());
			}*/
		    
			/*System.out.println("WINDOWS");
			System.out.println(windows.get(0).toString());
			System.out.println(windows.get(1).toString());
			for (int i = 0; i < windows.size(); i++)
			{
				System.out.println(i);
				System.out.println(windows.get(i).toString());
			}*/
			
			/*PrintWriter writer = new PrintWriter("windows.txt", "UTF-8");
			for (int i = 0; i < ancestral.size(); i++)
			{
				writer.println(ancestral.get(i).toString());
			}
			writer.close();*/
			
		    
		} 
		catch (IOException e) 
		{
			String msg = "VCF file not found";
			throw new FileParsingException(log, msg);
		}		
	}
	
	public boolean validAncestralData(String info)
	{
		boolean valid = false;
		if (info.contains("AA="))
		{
			List<String> a = new ArrayList<String>(Arrays.asList(info.split("AA=")));
			String allele = a.get(1).substring(0, 1).toUpperCase();
			if (allele.equals("A") || allele.equals("C") || allele.equals("G") || allele.equals("T") || allele.equals("?"))
				valid = true; 
		}
		return valid; 
	}
	
	public String ancestralAllele(List<String> line)
	{
		List<String> info = new ArrayList<String>(Arrays.asList(line.get(7).split("AA=")));
		String allele = info.get(1).substring(0, 1).toUpperCase();
		if (allele.equals("?"))
		{
			List<String> a = new ArrayList<String>(Arrays.asList(info.get(1).substring(2).split("\\|")));
			if (a.get(0).length() < a.get(1).length() || a.get(0) == "-") //if ancestral version is shorter
			{
				if (line.get(3).length() < line.get(4).length()) //return min(a0, a1)
					return line.get(3); 
				else
					return line.get(4); 
			}
			else //if ancestral version is longer
			{
				if (line.get(3).length() > line.get(4).length()) //return max(a0, a1)
					return line.get(3);
				else
					return line.get(4);
			}
		}
		else
			return allele;
	}
	
	public static void main(String[] args) 
	{
		VCFParser parser = new VCFParser("minime.vcf", 21, new Log());
		try {
			parser.parseVCF(1000000, true);
		} catch (FileParsingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public List<Window> getWindows()
	{
		return windows; 
	}
	public List<Window> getAncestral()
	{
		return ancestral;
	}
	public Individual[] getIndividuals()
	{
		return individuals; 
	}
	
}
