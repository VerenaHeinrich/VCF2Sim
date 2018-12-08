import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public class VCF{
	
	private String CHR = "CHR";
	private String POS = "POS";
	private String SNPID = "SNPID";
	private String REF = "REF";
	private String ALT = "ALT";
	private String QUAL = "QUAL";
	private String FILTER = "FILTER";
	private String NAMES = "NAMES";
	private String GENOTYPE = "GENOTYPE";
	private String FORMAT_GT = "FORMAT_GT";
	
	private String string_split = "-";
	
	public void get_individuum(	String[] list,
								String dir,
								HashMap<String,ArrayList<String>> info,
								HashMap<String, ArrayList<String>> genotypes,
								HashMap<String, String> reference_genotype,
								HashMap<String, ArrayList<ArrayList<Integer>>> region) throws IOException{
		
		for(String vcf:list){
			System.out.print("\nRead "+vcf+"\n");
						
			/**-------------------------------------**/
			/**Initialize attributes for vcf file**/
			
			HashMap<String, ArrayList<String>> current = new HashMap<String, ArrayList<String>>();
			ArrayList<String> names = new ArrayList<String>();
			ArrayList<Double> synonymous = new ArrayList<Double>();
			ArrayList<Double> nonsynonymous = new ArrayList<Double>();
			ArrayList<Double> transitions = new ArrayList<Double>();
			ArrayList<Double> transversions = new ArrayList<Double>();
			ArrayList<Double> heterozygous = new ArrayList<Double>();
			ArrayList<Double> homozygous = new ArrayList<Double>();
			ArrayList<Double> dbsnp = new ArrayList<Double>();
			ArrayList<Double> number_of_rows = new ArrayList<Double>();
			ArrayList<ArrayList<Double>> het_freq = new ArrayList<ArrayList<Double>>();
			string_split = "-";
			
			/**-------------------------------------**/
			/**Parse vcf file**/

			FileReader file = null;
			try {
				file = new FileReader(dir+vcf);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			BufferedReader in = new BufferedReader(file);

			int row_count = 0;
			String row = null;
			read: while((row = in.readLine()) != null){
				
				/**Count lines in vcf file:**/
				row_count++;
				if((row_count % 5000) == 0)	System.out.print("Processing row: "+row_count+"\n");
				
				/** header:**/
				if(row.startsWith("##"))	continue;
				if(row.startsWith("#CHROM")){
					String[] first_line = row.split("\t");
					for(int i = 9;i<first_line.length;i++) names.add(first_line[i]);
				}else{
					
					/**Parse Current VCF Line**/
					current = get_line(row, true);
					if(current.isEmpty()) continue;
					
					/** check if region is within consensus region **/
					if(check_region(region,  current.get(CHR).get(0),  Integer.valueOf(current.get(POS).get(0))) == false) continue;
					
					/**Skip INDELS:**/
					for(String e:current.get(REF)) if(	e.matches(".*[A-Za-z]{2,}.*") |
														e.equals("-") |
														e.equals("N") |
														e.equals(".") |
														e.contains("*")
													) continue read;
					for(String e:current.get(ALT)) if(	e.matches(".*[A-Za-z]{2,}.*") |
														e.equals("-") |
														e.equals("N")|
														e.equals(".")|
														e.contains("*")
													) continue read;
					
					/**--------------------------------------------**/
					/**Get additional Information for the vcf file:**/		
					if(current.get(ALT).size() == 1){
							
						/**Nonsynonymous/Synonymous Ratio**/
						get_synonymous(current,synonymous,nonsynonymous);
								
						/**Transitions/Transversions Ratio**/
						get_trans(current,transitions,transversions);
										
						/**Heterozygous/Homozygous Ratio**/
						get_heterozygous_homozygous(current, heterozygous, homozygous);
									
						/**DBSNP Ratio**/
						get_dbsnp_entries(current, dbsnp, number_of_rows);
									
						/**Heterozygous Allele Frequency**/
						get_heterozygous_allele_frequency(current, het_freq);
					}
							
					/**Skip chromosomes:**/
					if(	current.get(CHR).get(0).matches(".*[YXMMT].*") ) 	continue;
					
					
					/**Add Genotype Information:**/
					String key = current.get(CHR).get(0)+"\t"+current.get(POS).get(0);
					if(!genotypes.containsKey(key)){genotypes.put(key, new ArrayList<String>());
					
						/**Fill up Genotypes for individuals which were not present in earlier vcf file:**/
						for(int j=0; j<(info.get(NAMES).size());j++)	{
							genotypes.get(key).add(current.get(REF).get(0)+current.get(REF).get(0));
						}
					}
					
					for(int i = 0;i<current.get(FORMAT_GT).size();i++){
						genotypes.get(key).add(current.get(GENOTYPE).get(i));
					}
					
					/**Add reference genotype (== 0|0):**/
					if(!reference_genotype.containsKey(key)) reference_genotype.put(key, current.get(REF).get(0)+current.get(REF).get(0));
						
				}
			}
			
			/**-------------------------------------**/
			/**Collect individual names**/
			if(names.isEmpty()){
				for(int i = 0;i<current.get(FORMAT_GT).size();i++) info.get(NAMES).add("Sample"+(info.get(NAMES).size()+1));
			}else info.get(NAMES).addAll(names);
			
			/**-------------------------------------**/
			/**Add vcf attributes to hash**/
			for(int i = 0;i<current.get(FORMAT_GT).size();i++){
				/**-------------------------------------**/
				/**Add Nonsynonymous/Synonymous Ratio**/
				info.get("SYNONYMOUSRATIO").add(String.valueOf(nonsynonymous.get(i)/synonymous.get(i)));
			
				/**-------------------------------------**/
				/**Add Transitions/Transversins Ratio**/
				info.get("TRANSRATIO").add(String.valueOf(transitions.get(i)/transversions.get(i)));
				
				/**-------------------------------------**/
				/**Add Heterozygous/Homozygous Ratio**/
				info.get("HETHOMRATIO").add(String.valueOf(heterozygous.get(i)/homozygous.get(i)));
				
				/**-------------------------------------**/
				/**Add Percentage in dbSNP**/
				if(dbsnp.get(i) == 0.0) info.get("DBSNP").add("NaN");
				else info.get("DBSNP").add(String.valueOf(dbsnp.get(i)/number_of_rows.get(i)));
				
				/**-------------------------------------**/
				/**Add Variance of heterozygous allele frequency**/
				info.get("VARHET").add(String.valueOf(getVariance(het_freq.get(i))));
				
			}

		}
		System.out.print("("+info.get(NAMES).size()+" individual(s))\n");
	}
	
	/**-***************************************************************************************-**/
	/**Funktions:**/
	
	public void get_heterozygous_allele_frequency(HashMap<String, ArrayList<String>> current, ArrayList<ArrayList<Double>> het_freq){
		
		/**Initialize:**/
		if(het_freq.isEmpty()) for(@SuppressWarnings("unused") String e:current.get(FORMAT_GT)) het_freq.add(new ArrayList<Double>());
		
		int ref_ac = 0;		
		int alt_ac = 0;		
		int covered_individual_count = 0;
		
		if(current.get("SNPID").get(0).startsWith("rs")){
			for(int i = 0;i<current.get(FORMAT_GT).size();i++){
				
				/**Count covered Positions:**/
				if(!current.get(FORMAT_GT).get(i).matches(".*\\..*"))	covered_individual_count++;

				/**Just Calculate for Calls with dbSNP entry:**/			
				if(		current.get(FORMAT_GT).get(i).contains("1")|
						current.get(FORMAT_GT).get(i).contains("2")|
						current.get(FORMAT_GT).get(i).contains("3")
						)
				{		
					/**Just for heterozygous calls:**/
					if(		!(	current.get(FORMAT_GT).get(i).matches("1.1")|
								current.get(FORMAT_GT).get(i).matches("2.2")|
								current.get(FORMAT_GT).get(i).matches("3.3")
							)
					){						
						/**Take either the AD or the DP4 flag::**/
						if(current.containsKey("DP4")){
							ref_ac = Integer.valueOf(current.get("DP4").get(0))+Integer.valueOf(current.get("DP4").get(1));
							alt_ac = Integer.valueOf(current.get("DP4").get(2))+Integer.valueOf(current.get("DP4").get(3));
						}else if(current.containsKey("FORMAT_AD")){
							
							if(current.get("FORMAT_AD").get((covered_individual_count-1)*2).equals(".")) continue;
							if(current.get("FORMAT_AD").get((covered_individual_count-1)*2).equals(".")) continue;
							
							ref_ac = Integer.valueOf(current.get("FORMAT_AD").get((covered_individual_count-1)*2));
							alt_ac = Integer.valueOf(current.get("FORMAT_AD").get(((covered_individual_count-1)*2)+1));
						}
							
						if(ref_ac + alt_ac >= 20){
							
							/**downsample to 20x:**/
						   	char[] ref_array = new char[ref_ac];
						   	Arrays.fill(ref_array,'a');
						   	
						   	char[] alt_array = new char[alt_ac];
						   	Arrays.fill(alt_array,'b');
						    	
						    String ac = new String(ref_array)+new String(alt_array);
						    List<String> ref_list = Arrays.asList( ac.split("") );
						    Collections.shuffle(ref_list);
						        
						    String res = new String();
						    for(String r:ref_list){
						      	res = res+r;
						    }
						        
						    res = res.substring(0,20);
						    res = res.replace("b", "");
						        
						    int ref_length = res.length();		    	
						
					    	het_freq.get(i).add(Double.valueOf(ref_length)/20.0);
						}
					}
				}
			}
		}
	}
	
	public void get_dbsnp_entries(HashMap<String, ArrayList<String>> current, ArrayList<Double> dbsnp, ArrayList<Double> number_of_rows){
		
		/**Initialize:**/
		if(dbsnp.isEmpty()) for(@SuppressWarnings("unused") String e:current.get(FORMAT_GT)) dbsnp.add(0.0);
		if(number_of_rows.isEmpty()) for(@SuppressWarnings("unused") String e:current.get(FORMAT_GT)) number_of_rows.add(0.0);
		
		for(int i = 0;i<current.get(FORMAT_GT).size();i++){
					
			/**Just Calculate for Calls**/
			if(		current.get(FORMAT_GT).get(i).contains("1")|
					current.get(FORMAT_GT).get(i).contains("2")|
					current.get(FORMAT_GT).get(i).contains("3")
			){
				number_of_rows.set(i, number_of_rows.get(i)+1.0);
				
				if(current.get("SNPID").get(0).startsWith("rs")) dbsnp.set(i, dbsnp.get(i)+1.0);
			}
		}	
	}
	
	public void get_heterozygous_homozygous(HashMap<String, ArrayList<String>> current, ArrayList<Double> heterozygous,ArrayList<Double> homozygous ){
		
		/**Initialize:**/
		if(heterozygous.isEmpty()) for(@SuppressWarnings("unused") String e:current.get(FORMAT_GT)) heterozygous.add(0.0);
		if(homozygous.isEmpty()) for(@SuppressWarnings("unused") String e:current.get(FORMAT_GT)) homozygous.add(0.0);
				
		for(int i = 0;i<current.get(FORMAT_GT).size();i++){
			
			/**Just Calculate for Calls**/
			if(		current.get(FORMAT_GT).get(i).contains("1")|
					current.get(FORMAT_GT).get(i).contains("2")|
					current.get(FORMAT_GT).get(i).contains("3")
			){
				if(		current.get(FORMAT_GT).get(i).matches("1.1")|
						current.get(FORMAT_GT).get(i).matches("2.2")|
						current.get(FORMAT_GT).get(i).matches("3.3")
				)	homozygous.set(i, homozygous.get(i)+1.0);	
				else heterozygous.set(i, heterozygous.get(i)+1.0);	
			}
		}
	}

	public void get_trans(HashMap<String, ArrayList<String>> current,ArrayList<Double> transitions,ArrayList<Double> transversions){
		
		/**Initialize:**/
		if(transitions.isEmpty()) for(@SuppressWarnings("unused") String e:current.get(FORMAT_GT)) transitions.add(0.0);
		if(transversions.isEmpty()) for(@SuppressWarnings("unused") String e:current.get(FORMAT_GT)) transversions.add(0.0);
		
		for(int i = 0;i<current.get(FORMAT_GT).size();i++){
			
			/**Just Calculate for Calls**/
			if(		current.get(FORMAT_GT).get(i).contains("1")|
					current.get(FORMAT_GT).get(i).contains("2")|
					current.get(FORMAT_GT).get(i).contains("3")
			){
				if(current.get(REF).get(0).equals(current.get(ALT).get(0)))	transitions.set(i, transitions.get(i)+1.0);
				
				else if(current.get(REF).get(0).equals("A")){
					 if(current.get(ALT).get(0).equals("C"))		transversions.set(i, transversions.get(i)+1.0);
					 else if(current.get(ALT).get(0).equals("G"))	transitions.set(i, transitions.get(i)+1.0);
					 else if(current.get(ALT).get(0).equals("T"))	transversions.set(i, transversions.get(i)+1.0);
				}
				else if(current.get(REF).get(0).equals("C")){
					 if(current.get(ALT).get(0).equals("A"))		transversions.set(i, transversions.get(i)+1.0);
					 else if(current.get(ALT).get(0).equals("G"))	transversions.set(i, transversions.get(i)+1.0);
					 else if(current.get(ALT).get(0).equals("T"))	transitions.set(i, transitions.get(i)+1.0);
				}
				else if(current.get(REF).get(0).equals("G")){
					 if(current.get(ALT).get(0).equals("A"))		transitions.set(i, transitions.get(i)+1.0);
					 else if(current.get(ALT).get(0).equals("C"))	transversions.set(i, transversions.get(i)+1.0);
					 else if(current.get(ALT).get(0).equals("T"))	transversions.set(i, transversions.get(i)+1.0);
				}
				else if(current.get(REF).get(0).equals("T")){
					 if(current.get(ALT).get(0).equals("A"))		transversions.set(i, transversions.get(i)+1.0);
					 else if(current.get(ALT).get(0).equals("G"))	transversions.set(i, transversions.get(i)+1.0);
					 else if(current.get(ALT).get(0).equals("C"))	transitions.set(i, transitions.get(i)+1.0);
				}
			}
		}
	}
	
	public void get_synonymous(HashMap<String, ArrayList<String>> current,ArrayList<Double> synonymous,ArrayList<Double> nonsynonymous){
		
		/**Initialize:**/
		if(synonymous.isEmpty()) for(@SuppressWarnings("unused") String e:current.get(FORMAT_GT)) synonymous.add(0.0);
		if(nonsynonymous.isEmpty()) for(@SuppressWarnings("unused") String e:current.get(FORMAT_GT)) nonsynonymous.add(0.0);
		
		if(current.containsKey("EFFECT")){
			if(current.get("EFFECT").get(0).equals("synonymous")){
				for(int i = 0;i<current.get(FORMAT_GT).size();i++){
					
					/**Just Calculate for Calls**/
					if(		current.get(FORMAT_GT).get(i).contains("1")|
							current.get(FORMAT_GT).get(i).contains("2")|
							current.get(FORMAT_GT).get(i).contains("3")
					){
						synonymous.set(i, synonymous.get(i)+1.0);
					}
				}
			}			
			else if(	current.get("EFFECT").get(0).equals("stopgain") |
						current.get("EFFECT").get(0).equals("stoploss") | 
						current.get("EFFECT").get(0).equals("nonsynonymous")|
						current.get("EFFECT").get(0).equals("missense")
			){						
				for(int i = 0;i<current.get("FORMAT_GT").size();i++){
					
					/**Just Calculate for Calls**/
					if(		current.get(FORMAT_GT).get(i).contains("1")|
							current.get(FORMAT_GT).get(i).contains("2")|
							current.get(FORMAT_GT).get(i).contains("3")
					){
						nonsynonymous.set(i, nonsynonymous.get(i)+1.0);
					}
				}
			}
		}
	}
	public HashMap<String, ArrayList<String>> get_line(String row, boolean all) {
		
		HashMap<String, ArrayList<String>> line = new HashMap<String, ArrayList<String>>();
		if(row.startsWith("#"))	return line;	
		
		String[] current =  row.split("\t");  
				
		String[] info = null;
		String[] format = null;
		
		line.put(CHR, new ArrayList<String>());
		line.put(POS, new ArrayList<String>());
		line.put(REF, new ArrayList<String>());
		line.put(ALT, new ArrayList<String>());
		line.put(GENOTYPE, new ArrayList<String>());
		
		/**-------------------------------------**/
		/**Define fixed Attributes of a VCF file**/
	
		line.get(CHR).add(current[0].replaceAll("chr|Chr", ""));
		line.get(POS).add(current[1]);
		for(String e: current[3].split(",")) line.get(REF).add(e.toUpperCase()); 
		for(String e: current[4].split(",")) line.get(ALT).add(e.toUpperCase()); 
		
		/**-------------------------------------**/
		/**Define variable Attributes of a VCF file**/
					
		if(all == true){	
			line.put(SNPID, new ArrayList<String>());
			line.put(QUAL, new ArrayList<String>());
			line.put(FILTER, new ArrayList<String>());
			
			line.get(SNPID).add(current[2]);
			line.get(QUAL).add(current[5]);
			line.get(FILTER).add(current[6]);
			
			/**INFO column:**/
			if(current[7].contains(";")){
				
				info = current[7].split(";|=");
				for(int i = 0;i<(info.length-1);i=i+2){
					line.put(info[i], new ArrayList<String>());
					
					for(String a:info[i+1].split(",")) {
						line.get(info[i]).add(a);
					}
				}
			}else line.put(current[7], new ArrayList<String>());
		}
		
		/**FORMAT column:**/
		if(current[8].contains(":")){
			format = current[8].split(":");
			for(String e: format)	line.put(("FORMAT_"+e), new ArrayList<String>());
					
			/**Individual column:**/
			for(int i = 9;i<current.length;i++){
				String[] ind = current[i].split(":");
				
				for(int j = 0;j<ind.length;j++) {
					if(format[j].equals("AD")){
						if(ind[j].equals(".")){
							ind[j] = "0,0";
						}
					}
					for(String a:ind[j].split(",")) {
						line.get("FORMAT_"+format[j]).add(a);
						
					}
				}
			}
		}else {
			line.put(("FORMAT_GT"), new ArrayList<String>());
			
			/**Individual column:**/
			for(int i = 9;i<current.length;i++){
				String[] ind = current[i].split(":");
				line.get("FORMAT_GT").add(ind[0]);
			}
		}
		
		/**Save time by remembering old pattern:**/
		HashMap<String, String> already_seen_genotype = new HashMap<String, String>();
		
		/**Get Actual Genotype:**/
		if(string_split.equals("-")){
			if(line.get(FORMAT_GT).get(0).contains("|")) string_split = "\\|";
			else if(line.get(FORMAT_GT).get(0).contains("/")) string_split = "\\/";
			else string_split = "";
		}
		
		/**Get Actual Genotype:**/
		for(String a:line.get(FORMAT_GT)){
			
			/**Get string split pattern: **/
			if(a.matches(".*\\|.*")) string_split = "\\|";
			else if (a.matches(".*\\/.*"))  string_split = "\\/";
			else{
				
				/**In hemizygote chrX stat:**/
				a = a+"|"+a;
				string_split = "\\|";
			}
			
			/**Just put in already known genotype:**/
			if(already_seen_genotype.containsKey(a))	{
				line.get(GENOTYPE).add(already_seen_genotype.get(a));
				continue;
			}

			List<String> tmp = new ArrayList<String>();
			for(String e:a.split(string_split)){
				if(e.equals("0")) tmp.add(line.get(REF).get(0));
				else if(e.contains(".")) tmp.add(line.get(REF).get(0));//tmp.add(".");
				else tmp.add(line.get(ALT).get(Integer.valueOf(e)-1));
			}

			Collections.sort(tmp);
			already_seen_genotype.put(a, tmp.get(0)+tmp.get(1));
			line.get(GENOTYPE).add(already_seen_genotype.get(a));			
		}
		return line;
	}
	
	/** check if positon is within defined consensus region **/
	Boolean check_region(HashMap<String, ArrayList<ArrayList<Integer>>> region, String chr, Integer pos){
		
		Boolean check=false;
		// go through all regions in region:
		for(ArrayList<Integer> this_region:region.get(chr)){
			if (this_region.get(0) <= pos & this_region.get(1) >= pos){
				check = true;
				break;
			}

		}
		return check;
	}
	
	void update_genotypes(	HashMap<String, ArrayList<String>> genotypes,
							HashMap<String,ArrayList<String>> info,
							HashMap<String, String> reference_genotype){
		
		int number_of_individuals = info.get("NAMES").size();
		for(String id:genotypes.keySet()){
			if(genotypes.get(id).size() != number_of_individuals){
				for(int i = genotypes.get(id).size(); i<number_of_individuals;i++){
					genotypes.get(id).add(reference_genotype.get(id));
				}
			}
		}
	}
	
	double getMean(ArrayList<Double> data)
    {
        double sum = 0.0;
        double size = data.size();
        
        for(double a : data)
            sum += a;
            return sum/size;
    }
	
    double getVariance(ArrayList<Double> data)
    {
        ArrayList<Double> data_sqrt = new ArrayList<Double>();
        
        for(double a:data){
        	data_sqrt.add(a*a);
        }
        
        return (getMean(data_sqrt)-(getMean(data)*getMean(data))); 		
    }
	/**-***************************************************************************************-**/
}