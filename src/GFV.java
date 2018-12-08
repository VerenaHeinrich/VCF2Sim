import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class GFV{
	
	private ArrayList<String> GENOTYPES = new ArrayList<String>();
	public static String[] line;

	public Double get_gfv(String gfv_file, HashMap<String, HashMap<String, Double>> vector, double max_value) throws IOException{
		
		System.out.print("Read genotype frequency vector\n");
		
		/**--------------------------**/
		/**All 10 Possible Genotypes:**/
		GENOTYPES.add("AA");
		GENOTYPES.add("AC");
		GENOTYPES.add("AG");
		GENOTYPES.add("AT");
		GENOTYPES.add("CC");
		GENOTYPES.add("CG");
		GENOTYPES.add("CT");
		GENOTYPES.add("GG");
		GENOTYPES.add("GT");
		GENOTYPES.add("TT");
	
		/**Initializations:**/
		String[] current = null;
		String key = null;
		String genotype_key = "";
		double tmp_max_value = 0.0;
		double value_from_gfv = 0;

		/**--------------**/
		/**Get max value:**/
		FileReader file_tmp = null;
		try {
			file_tmp = new FileReader(gfv_file);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		BufferedReader in_tmp = new BufferedReader(file_tmp);
		
		String row_tmp = null;
		while((row_tmp = in_tmp.readLine()) != null) {
			
			/**Current row:**/
			current =  row_tmp.split("\t");  
			key = current[0]+"\t"+current[1];
						
			/**get maximum individual value from reference set:**/ 
			tmp_max_value = 0.0;
			for(int i=0;i<10;i++){
				value_from_gfv = Double.valueOf(current[4+i]);
				tmp_max_value += value_from_gfv;
			}
			if(max_value <= tmp_max_value) max_value = tmp_max_value;
		}
		
		/**--------------**/
		/**Read GFV file:**/
		FileReader file = null;
		try {
			file = new FileReader(gfv_file);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		BufferedReader in = new BufferedReader(file);
		
		String row = null;
		while((row = in.readLine()) != null) {
			
			/**Current row:**/
			current =  row.split("\t");  
			key = current[0]+"\t"+current[1];
			
			HashMap<String,Double> init_gfv = new HashMap<String, Double>();

			/**Assign frequencies at each position:**/
			for(int i=0;i<10;i++)	{
				value_from_gfv = Double.valueOf(current[4+i]);
				genotype_key = GENOTYPES.get(i);
				
				if(value_from_gfv > 0.0)	{
					init_gfv.put(genotype_key, value_from_gfv/max_value);
				}else{
					init_gfv.put(genotype_key, 0.0);
				}
			}
			vector.put(key, init_gfv);
		}
		return(max_value);
	}
	
	public void update_gfv(HashMap<String, ArrayList<String>> genotypes, HashMap<String, HashMap<String, Double>> vector, HashMap<String, String> reference_genotype, double max_value){
		
		System.out.print("\nUpdate genotype frequency vector: \n");
		
		/**Initialize variable to count every position, which is not covered in the GFV file:**/
		int new_positions = 0;
		
		for(String e:genotypes.keySet()){
			if(!vector.containsKey(e)) {
				
				/**Sum up new positions:**/
				new_positions++;
				
				/**Reference Genotye:**/
				String ref_gt = reference_genotype.get(e);
				
				/**Initialize gfv vector:**/
				vector.put(e, new HashMap<String,Double>());
				
				/**Initialize frequencies at each position:**/
				for(int i=0;i<10;i++){
					if(GENOTYPES.get(i).equals(ref_gt)){
						vector.get(e).put(GENOTYPES.get(i), max_value);
					}else{
						vector.get(e).put(GENOTYPES.get(i), 0.0);
					}
				}
				
				/**Set genotype from sample to 1.0:**/
				for(int i=0;i<10;i++){
					for(int j=0;j<genotypes.get(e).size();j++){
						if(GENOTYPES.get(i).equals(genotypes.get(e).get(j)) & !GENOTYPES.get(i).equals(ref_gt)){
							vector.get(e).put(GENOTYPES.get(i),	1.0);	//just initialize this genotype with 1 individual
						}
					}
				}
				
				/**Normalization:**/
				for(int i=0;i<10;i++)	vector.get(e).put(GENOTYPES.get(i), vector.get(e).get(GENOTYPES.get(i))/(1.0+max_value));
								
			}
		}		
		System.out.print(new_positions+" new positions in the sample(s)\n\n");
	}
	
	public HashMap<String, ArrayList<ArrayList<Integer>>> read_consensus(String consensus_file) throws IOException, InterruptedException{
		
		/**HashMap that lists the regions in the 1KGP consensus region file:**/
		HashMap<String, ArrayList<ArrayList<Integer>>> region = new HashMap<String, ArrayList<ArrayList<Integer>>>();
				
		/**Open genotype count file:**/
		FileReader file = null;
		try {
			file = new FileReader(consensus_file);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		BufferedReader in = new BufferedReader(file);

		String row = null;
		while((row = in.readLine()) != null){
			
			/**Split current line:**/
			line = row.split("\t");
			line[0] = line[0].replaceAll("chr|Chr", "");
			
			/**Initialize Reference Genotype:**/
			ArrayList<Integer> init = new ArrayList<Integer>();
			init.add(Integer.valueOf(line[1]));
			init.add(Integer.valueOf(line[2]));
				
			if(!region.containsKey(line[0])){
				region.put(line[0], new ArrayList<ArrayList<Integer>>());
				region.get(line[0]).add(0, init);
			}else{
				region.get(line[0]).add(init);
			}
			
		}
		in.close();
		return region;
	}
}