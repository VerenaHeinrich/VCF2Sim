import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class main{
	
	/**-***************************************************************************************-**/

	/**-----------------------------------------------**/
	/**Initialize additional information for vcf file:**/
	private static String NAMES = "NAMES";
	private static String SYNONYMOUSRATIO = "SYNONYMOUSRATIO";
	private static String TRANSRATIO = "TRANSRATIO";
	private static String HETHOMRATIO = "HETHOMRATIO";
	private static String DBSNP = "DBSNP";
	private static String VARHET = "VARHET";

	/**-***************************************************************************************-**/

	public static void main(String[] args) throws Exception{
		
		/**---------------**/
		/**Start Time Count:**/
		long startTime = System.nanoTime();
		
		/**Command Line Options**/
		Commands com = new Commands();
		ArrayList<String> commands = com.get_commands(args);
		String dir = commands.get(0);
		String pop = commands.get(1);
		String out = commands.get(2);
		String metric = commands.get(3);
		String gfv_file = commands.get(4);

		/**Get .vcf files**/
		String[] vcf_list;
		if(dir.endsWith(".vcf") |dir.endsWith(".vcf.txt")){
			String[] tmp = dir.split("/");
			vcf_list = tmp[(tmp.length-1)].split("\n");
			dir = dir.replaceAll(vcf_list[0],"");
		}else{
			vcf_list = get_file_names(dir);
		}

		/**INFO for vcf files:**/
		HashMap<String, ArrayList<String>> info = new HashMap<String, ArrayList<String>>();
		
		/**Attributes of INFO-hash:**/
		info.put(NAMES, new ArrayList<String>());
		info.put(SYNONYMOUSRATIO, new ArrayList<String>());
		info.put(TRANSRATIO, new ArrayList<String>());
		info.put(HETHOMRATIO, new ArrayList<String>());
		info.put(DBSNP, new ArrayList<String>());
		info.put(VARHET, new ArrayList<String>());
		
		/**-------------------------------------------------------------------------------**/
		/**Read file for Genotype Frequency Vector (GFV):**/
		double max_value = 0.0;
		
		HashMap<String, HashMap<String, Double>> vector = new HashMap<String, HashMap<String,Double>>();
		GFV gfv = new GFV();
		max_value = gfv.get_gfv(gfv_file,vector,max_value);

		/**-------------------------------------------------------------------------------**/
		/**Initialize Hash with Reference Genotypes:**/
		HashMap<String, String> reference_genotype = new HashMap<String, String>(); 
		
		/**-------------------------------------------------------------------------------**/
		/**read consensus region:**/
		HashMap<String, ArrayList<ArrayList<Integer>>> region = gfv.read_consensus("data/20110225.exome.consensus.bed");

		/**-------------------------------------------------------------------------------**/
		/**Parse vcf file**/
		HashMap<String, ArrayList<String>> genotypes = new HashMap<String, ArrayList<String>>();
		VCF parse = new VCF();			
		parse.get_individuum(vcf_list, dir, info, genotypes, reference_genotype, region);
		
		/**-------------------------------------------------------------------------------**/
		/**Update genotypes: (positions which are not covered in some individuals)**/
		parse.update_genotypes(genotypes, info, reference_genotype);
		
		/**-------------------------------------------------------------------------------**/
		/**update genotype frequency vector at positions which are not present in the gfv:**/
		gfv.update_gfv(genotypes, vector, reference_genotype, max_value);
		
		/**-------------------------------------------------------------------------------**/
		/**get best fitting background population from 1KGP:**/
		Compare compare = new Compare();
		if(pop == null){
			System.out.print("Check best fitting 1KG Background Population:\n------------------------------------------------------\n");
			ArrayList<String> names = get_population_names("data/first_individuals_all_pop_collapsed.vcf");
			
			@SuppressWarnings("unchecked")
			ArrayList<ArrayList<Double>> check_matrix = compare.compare_to_background(	(HashMap<String, ArrayList<String>>) genotypes.clone(),
																						"first_individuals_all_pop",
																						info,
																						vector, 
																						reference_genotype,
																						dir+out,
																						false,
																						metric
																					);
			pop=get_best_population(check_matrix, names, info);
		}
		
		/**-------------------------------------------------------------------------------**/
		System.out.print("Compare to Background Population: "+pop+"\n");

		/**-------------------------------------------------------------------------------**/
		/**Parse Background Population vcf file and compare with individuals:**/		
		@SuppressWarnings("unused")
		ArrayList<ArrayList<Double>> matrix = compare.compare_to_background(	genotypes,
																				pop,
																				info,
																				vector,
																				reference_genotype,
																				dir+out,
																				true,
																				metric
																			);
		
		/**---------------**/
		/**Print Out Time:**/
		long endTime = System.nanoTime();
		System.out.println("\n------------------------------------------------------\nFinished!\nRuntime: "+(endTime - startTime)/1000000000 + " seconds"); 
	}
	
	/**-***************************************************************************************-**/
	/**Funktions:**/
	
	/**-----------------------------**/
	/**Get all vcf files in a directory
	 * @throws Exception **/
	public static String[] get_file_names(String dir) throws Exception{
				
		File file = new File(dir);
		if(!file.isDirectory()) throw new Exception("Directory "+dir+"does not exist.\n");
		
		String[] list = file.list(new FilenameFilter() {
			public boolean accept(File d, String name) {	return (name.endsWith(".vcf")|name.endsWith(".vcf.txt"));	}
		});	
		
		java.util.Arrays.sort(list);
		
		return list;
	}
	
	/**----------------------------**/
	/**Get best fitting background population from similarity matrix:**/
	public static String get_best_population(	ArrayList<ArrayList<Double>> matrix,
												ArrayList<String> names,
												HashMap<String, ArrayList<String>> info
											){
		int[] index = new int[matrix.get(0).size()];
		int[] count = new int[(matrix.size()-matrix.get(0).size())];
		
		for(int i=0;i<(matrix.get(0).size());i++){
			
			/**Initialization:**/
			double tmp = 2;
			
			for(int j=0;j<(matrix.size()-matrix.get(0).size());j++){
				if(matrix.get(j).get(i) <= tmp) {
					tmp = matrix.get(j).get(i);
					index[i] = j; 
				}
			}
			count[index[i]]++;
			System.out.print("Best fitting background population for "+info.get(NAMES).get(i)+" is: "+names.get(index[i])+"\n");
		}
		System.out.print("\n");
		
		/**get best fitting background population:**/
		int tmp = 0;
		for(int j=0;j<(matrix.size()-matrix.get(0).size());j++){
			if(count[j]>tmp){
				tmp = j;
			}
		}
		return names.get(tmp);
	}
	public static ArrayList<String> get_population_names(String input) throws IOException{
		ArrayList<String> names = new ArrayList<String>();
		
		FileReader file = null;
		try {
			file = new FileReader(input);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		BufferedReader in = new BufferedReader(file);

		String row = null;
		while((row = in.readLine()) != null){
			
			if(row.startsWith("#CHROM")){
				String[] first_line = row.split("\t");
				for(int i = 9;i<first_line.length;i++) names.add(first_line[i]);
			}
			
		}
		in.close();
		return names;
	}
	
	/**-***************************************************************************************-**/
}

