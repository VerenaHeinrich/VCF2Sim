import java.util.ArrayList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;


public class Commands{
	public static String dir = null;
	public static String out = null;
	public static String pop = null;
	public static String metric = "rare";
	public static String gfv_file = "data/GFV_consensus.txt";

	public ArrayList<String> get_commands(String[] args) throws Exception{
		try {
		      Options opt = new Options();
		      opt.addOption("d","dir", true, "Directory to vcf files");
		      opt.addOption("p","pop",true,"1000Genomes background population (if not stated, the most likely will be selected.)");
		      opt.addOption("o","out",true,"Name of output file");
		      opt.addOption("m","metric",true,"Metric (Sensitive for 'rare' or 'common' variants. You can also use a normal 'hamming' distance. default = 'rare')");
		      opt.addOption("g","gfv",true,"Name of optional GFV file. default build on 1KGP ('/VCF2Sim/data/GFV_consensus.txt')");

		      GnuParser parser = new GnuParser();
		      CommandLine cl = parser.parse(opt, args);
		      
		      if(cl.hasOption('d')){dir = cl.getOptionValue('d');}
		      else{
			         HelpFormatter f = new HelpFormatter();
			         f.printHelp("VCF2Sim.jar", opt);
			         System.exit(0);
			      }	
		      if(cl.hasOption('o')){out = cl.getOptionValue('o');}
		      else{
			         HelpFormatter f = new HelpFormatter();
			         f.printHelp("VCF2Sim.jar", opt);
			         System.exit(0);
			      }	
		      
		      /****/
		      if(cl.hasOption('p')){pop = cl.getOptionValue('p');}
		      if(cl.hasOption('m')){metric = cl.getOptionValue('m');}
		      if(cl.hasOption('g')){gfv_file = cl.getOptionValue('g');}

		      if(!metric.equals("rare") & !metric.equals("common") & !metric.equals("hamming")){
		    	  		    	  
		    	  HelpFormatter f = new HelpFormatter();
			      f.printHelp("VCF2Sim.jar", opt);
			      
		    	  throw new Exception("Metric "+metric+" is not an option.\n");
		      }
		      
		}catch (ParseException e) {
		      e.printStackTrace();
		}
		
		ArrayList<String> ret = new ArrayList<String>();
		ret.add(dir);
		ret.add(pop);
		ret.add(out);
		ret.add(metric);
		ret.add(gfv_file);
		
		return ret;
	}
}