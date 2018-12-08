import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public class Compare{
	
	/**-***************************************************************************************-**/
	
	/**-----------------------------------**/
	/**Initialize Attributes for vcf line:**/
	private String CHR = "CHR";
	private String POS = "POS";
	private String REF = "REF";
	private String NAMES = "NAMES";
	private String GENOTYPE = "GENOTYPE";
	
	/**-----------------------------------**/
	/**Initialize current line in vcf file**/
	private HashMap<String, ArrayList<String>> current = new HashMap<String, ArrayList<String>>();
	
	/**-------------------------------------------------**/
	/**Initialize values from genotype frequency vector:**/
	private Double value_from_sample = 0.0;
	private Double value_from_reference = 0.0;
	private Double value_norm = 0.0;
	
	/**---------------------**/
	/**Initialize Genotypes:**/
	private ArrayList<String> current_genotype = new ArrayList<String>();
	private String ref_genotype = "";
	private String genomes_gt = "";
	private String sample_gt = "";
	
	/**--------------------------------**/
	/**Initialize ID = 'chr\tposition':**/
	private String id = "";
	private HashMap<String, Double> vector_at_id = new HashMap<String, Double>();
	private ArrayList<String> genotypes_at_id = new ArrayList<String>();
	
	/**---------------------------------**/
	/**Initialize number of individuals:**/
	private int number_of_reference_individuals = 0;
	private int number_of_sample_individuals = 0;
	private int number_of_all_individuals = 0;
	
	/**------------------------------------------------------**/
	/**Initialize numerator and denumerator similarity matrix**/
	private ArrayList<ArrayList<Double>> matrix_numerator = new ArrayList<ArrayList<Double>>();
	private ArrayList<ArrayList<Double>> matrix_denumerator = new ArrayList<ArrayList<Double>>();
	private ArrayList<ArrayList<Double>> matrix = new ArrayList<ArrayList<Double>>();

	/**---------------------------------------**/
	/**Initialize Hash for already seen values**/
	private HashMap<String, Double> already_seen = new HashMap<String, Double>();

	/**-***************************************************************************************-**/

	public ArrayList<ArrayList<Double>> compare_to_background(	HashMap<String, ArrayList<String>> genotypes,
																String pop,
																HashMap<String, ArrayList<String>> info,
																HashMap<String, HashMap<String, Double>> vector,
																HashMap<String, String> reference_genotype,
																String out,
																boolean print,
																String metric
															) throws Exception{
		
		/**-------------------------------------**/
		/**clear matrices**/
		matrix_numerator.clear();
		matrix_denumerator.clear();
		matrix.clear();
		
		/**-------------------------------------**/
		/**Parse vcf file**/

		VCF parse = new VCF();			
		FileReader file = null;
		try {
			file = new FileReader("data/"+pop+"_collapsed.vcf");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		BufferedReader in = new BufferedReader(file);
		String row = null;
		while((row = in.readLine()) != null){
			
			/**Parse Current VCF Line**/
			current = parse.get_line(row, false);
			if(current.isEmpty()) continue;
			
			/**---------------**/
			/**Initialization:**/
			id = current.get(CHR).get(0)+"\t"+current.get(POS).get(0);
			current_genotype = current.get(GENOTYPE);
			
			if(!vector.containsKey(id)) continue;

			vector_at_id = vector.get(id);
			
			ref_genotype = current.get(REF).get(0)+current.get(REF).get(0);	
			if(matrix.isEmpty()) init(info);
			
			/**-----------------------------**/
			/**Check if ID exist in samples:**/
			if(!genotypes.containsKey(id)) {	
				ArrayList<String> reference_list = new ArrayList<String>();
				
				/**Fill up Genotypes:**/
				for(@SuppressWarnings("unused") String s:info.get(NAMES))	reference_list.add(ref_genotype);
				genotypes_at_id = reference_list;		
						
			} else {
				genotypes_at_id = genotypes.get(id);
			}

			/**--------------------------------------------------------**/
			/**Go through every line and compare with each individuals:**/
			for(int i=0;i<number_of_reference_individuals;i++){	
				genomes_gt = current_genotype.get(i);	
				compare_innner_loop(genomes_gt, vector_at_id, i+0, genomes_gt.equals(ref_genotype),0,metric,id,1);
			}
			
			/**-----------------------------------------------------------------------------------**/
			/**Compare samples to samples:**/
			for(int i=0;i<number_of_sample_individuals;i++){
				sample_gt = genotypes_at_id.get(i);
				compare_innner_loop(sample_gt, vector_at_id, i+number_of_reference_individuals, sample_gt.equals(ref_genotype),(i+1),metric,id,2);
			}
			
			/**Remove already compared keys from the sample set:**/
			genotypes.remove(id);
		}
		in.close();
		
		/**-------------------------------------------------------------------------**/
		/**Go through every remaining line that is not covered by the reference set:**/
		for(String id:genotypes.keySet()){			
			vector_at_id = vector.get(id);
			genotypes_at_id = genotypes.get(id);
			ref_genotype = reference_genotype.get(id);
			value_from_reference = vector_at_id.get(ref_genotype);

			List<Double> vector_sort = new ArrayList<Double>();
			for(String e:vector_at_id.keySet())	{
				if(vector_at_id.get(e) > 0.0) vector_sort.add(vector_at_id.get(e));
			}
			Collections.sort(vector_sort);
					
			for(int i=0;i<number_of_sample_individuals;i++){				
				sample_gt = genotypes_at_id.get(i);
				value_from_sample = vector_at_id.get(sample_gt);
				
				if(metric.equals("rare"))value_norm = 2/(value_from_sample+value_from_reference);
//				else value_norm = value_from_sample+value_from_reference;
				else if(metric.equals("common")){
					value_norm = 0.0;
					for(int e=0;e<vector_sort.size();e++)	value_norm += vector_sort.get(e)*Math.log(vector_sort.get(e));
				}else if(metric.equals("hamming")){
					value_norm = 1.0;
				}
						
				for(int k=0;k<(number_of_reference_individuals);k++){	
					if(sample_gt.equals(ref_genotype)) continue;
					matrix_denumerator.get(k).set(i, matrix_denumerator.get(k).get(i) + value_norm);
				}
				compare_innner_loop(sample_gt, vector_at_id, i+number_of_reference_individuals, sample_gt.equals(ref_genotype),(i+1),metric,id,3);
			}
			
//			/**-----------------------------------------------------------------------------------**/
//			/**Compare samples to samples:**/
//			for(int i=0;i<number_of_sample_individuals;i++){
//				sample_gt = genotypes_at_id.get(i);
//				compare_innner_loop(sample_gt, vector_at_id, i+number_of_reference_individuals, sample_gt.equals(ref_genotype),(i+1),metric,id,4);
//			}
		}
		
		/**Print Out Results:**/
		calculate_similarity_matrix(matrix_numerator, matrix_denumerator, matrix);
		if(print){
			print_out(matrix, out, pop, info);
			print_out(matrix_numerator, out+"_equal", pop, info);
			print_out(matrix_denumerator, out+"_norm", pop, info);
		}
		return matrix;
	}
	
	/**-***************************************************************************************-**/
	/**Funktions:**/
	
	void compare_innner_loop(	String sample_gt_first, 
								HashMap<String, Double>  vector_at_id, 
								int index,
								boolean first_equals_ref, 
								int start,
								String metric,String id,
								int c){
		
		/**----------------------**/
		/**Inner loop comparison:**/
		
		/**value from genotype frequency vector for first genotype:**/
		double value_from_sample_1st = vector_at_id.get(sample_gt_first);
		double value_from_sample_2nd = 0.0;
		
		List<Double> vector_sort = new ArrayList<Double>();
		for(String e:vector_at_id.keySet())	{
			if(vector_at_id.get(e) > 0.0) vector_sort.add(vector_at_id.get(e));
		}
		Collections.sort(vector_sort);
				
		already_seen.clear();
		
		for(int n=start;n<number_of_sample_individuals;n++){	
			String sample_gt_2nd = genotypes_at_id.get(n);
			if(first_equals_ref & sample_gt_2nd.equals(ref_genotype))	continue;
			
			/**Just put in already known value:**/
			if(already_seen.containsKey(sample_gt_2nd)){
				matrix_denumerator.get(index).set(n, matrix_denumerator.get(index).get(n) + already_seen.get(sample_gt_2nd));
				if(sample_gt_first.equals(sample_gt_2nd)) {
					matrix_numerator.get(index).set(n, matrix_numerator.get(index).get(n) + already_seen.get(sample_gt_2nd));
				}
				continue;
			}
			
			value_from_sample_2nd = vector_at_id.get(sample_gt_2nd);
			
			if(metric.equals("rare")) value_norm = 2/(value_from_sample_1st+value_from_sample_2nd);
//			else value_norm = value_from_sample_1st+value_from_sample_2nd;
			else if(metric.equals("common")){
				value_norm = 0.0;
				for(int e=0;e<vector_sort.size();e++)	value_norm += vector_sort.get(e)*Math.log(vector_sort.get(e));
			}else if(metric.equals("hamming")){
				value_norm = 1.0;
			}
			
			if(Double.isInfinite(value_norm)) continue;

			already_seen.put(sample_gt_2nd, value_norm);
			matrix_denumerator.get(index).set(n, matrix_denumerator.get(index).get(n) + value_norm);
			
			if(sample_gt_first.equals(sample_gt_2nd)) {
				matrix_numerator.get(index).set(n, matrix_numerator.get(index).get(n) + value_norm);
			}					
		}
	}
	void init(HashMap<String, ArrayList<String>> info){
		
		/**---------------------------------**/
		/**Initialize number of individuals:**/
		
		number_of_reference_individuals = current.get(GENOTYPE).size();
		number_of_sample_individuals = info.get(NAMES).size();
		number_of_all_individuals = number_of_reference_individuals+number_of_sample_individuals;
		
		matrix_numerator = init_matrix(number_of_all_individuals,number_of_sample_individuals);
		matrix_denumerator = init_matrix(number_of_all_individuals,number_of_sample_individuals);
		matrix = init_matrix(number_of_all_individuals,number_of_sample_individuals);
	}
	ArrayList<ArrayList<Double>> init_matrix(int I,int J){
		ArrayList<ArrayList<Double>> matrix = new ArrayList<ArrayList<Double>>();
		
		for(int i=0;i<I;i++){
			matrix.add(new ArrayList<Double>());

			for(int j=0;j<J;j++) matrix.get(i).add(0.0);
		}
		return matrix;
	}
	void calculate_similarity_matrix (ArrayList<ArrayList<Double>> matrix_numerator, ArrayList<ArrayList<Double>> matrix_denumerator, ArrayList<ArrayList<Double>> matrix) throws IOException{
		
		/**----------------------------**/
		/**Calculate Similarity Matrix:**/
		
		for(int i=0;i<number_of_all_individuals;i++){
			for(int j=0;j<number_of_sample_individuals;j++){
				
				/**Similarity Matrix:**/
				if((i-number_of_reference_individuals)==j) {
					matrix.get(i).set(j,0.0);
					matrix_numerator.get(i).set(j,0.0);
					matrix_denumerator.get(i).set(j,0.0);
				}
				else if((i-number_of_reference_individuals)>j)  {
					matrix.get(i).set(j,matrix.get(j+number_of_reference_individuals).get(i-number_of_reference_individuals));
					matrix_numerator.get(i).set(j,matrix_numerator.get(j+number_of_reference_individuals).get(i-number_of_reference_individuals));
					matrix_denumerator.get(i).set(j,matrix_denumerator.get(j+number_of_reference_individuals).get(i-number_of_reference_individuals));
				}
				else matrix.get(i).set(j, 1-matrix_numerator.get(i).get(j)/ matrix_denumerator.get(i).get(j));
				
			}
		}
	}
	/**----------------------------**/
	/**Print Out similarity matrix:**/
	public void print_out(ArrayList<ArrayList<Double>> matrix, String out, String pop, HashMap<String, ArrayList<String>> info) throws IOException{
		System.out.print(out+"\n");
		File f = new File(out);
		FileWriter w = new FileWriter(f);
		StringBuffer o = new StringBuffer();
		
		o.append("##pop\tsample name:ti/tv-ratio:Percentage in dbsnp138:var(het):nonsyn/syn-ratio:het/hom\n");
		o.append("#"+pop+"\t");
		
		/**print out additional information:**/
		for(int i=0;i<matrix.get(0).size();i++){
			o.append(info.get(NAMES).get(i).toString());
			o.append(":"+info.get("TRANSRATIO").get(i).toString());
			o.append(":"+info.get("DBSNP").get(i).toString());
			o.append(":"+info.get("VARHET").get(i).toString());
			o.append(":"+info.get("SYNONYMOUSRATIO").get(i).toString());
			o.append(":"+info.get("HETHOMRATIO").get(i).toString());

			if(i<(matrix.get(0).size()-1)) o.append("\t");
			else o.append("\n");
		}
		
		/**print out NaN for not filled values:**/
		for(int i=0;i<(matrix.size()-info.get(NAMES).size());i++){
			for(int j=0;j<matrix.size();j++){
				o.append("NaN"+"\t");
			}
			o.append("\n");
		}
		
		/**print out similarity values:**/
		for(int i=0;i<matrix.get(0).size();i++){
			for(int j=0;j<matrix.size();j++){
				o.append(matrix.get(j).get(i).toString()+"\t");
			}
			o.append("\n");
		}
		
		w.write(o.toString());
		w.close();
		
	}
	
	/**-***************************************************************************************-**/
}