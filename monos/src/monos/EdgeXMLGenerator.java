package monos;

import java.io.*;
import java.util.*;

import beast.app.beauti.BeautiDoc;
import beast.app.util.Application;
import beast.app.util.OutFile;
import beast.app.util.XMLFile;
import beast.core.Description;
import beast.core.Input;
import beast.core.util.Log;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.Exponential;
import beast.math.distributions.LogNormalDistributionModel;
import beast.math.distributions.MRCAPrior;
import beast.math.distributions.Normal;
import beast.math.distributions.ParametricDistribution;
import beast.math.distributions.Uniform;
import beast.util.NexusParser;
import beast.util.TreeParser;

@Description("For getting information from the  GlottoLog newick tree file and produce edge.xml")
public class EdgeXMLGenerator extends beast.core.Runnable {
	final public Input<File> treeFileInput = new Input<>("glottologtree", "file containins glottolog newick tree file", new File("data/tree_glottolog_newick.txt"));
	final public Input<File> familyExclusionsInput = new Input<>("familyExclusions", "text file with GlottoLog codes of families that should not be included in the analysis (e.g., sign languages)", new File("data/familyExclusions.csv"));
	final public Input<File> languageExclusionsInput = new Input<>("languageExclusions", "text file with GlottoLog codes of languages that should not be included in the analysis", new File("data/languageExclusions.dat"));
	final public Input<File> cladesInput = new Input<>("cladesDir", "specify directory with clades over families, e.g., Americas", new File("data/clades"));
	final public Input<File> calibrationsInput = new Input<>("calibrations", "NEXUS file with time calibrations for MRCAs of clades", new File("data/calibrations.nex"));
	final public Input<File> originateCalibrationsInput = new Input<>("originateCalibrations", "NEXUS file with time calibrations for orignate of MRCA of clades (including singletons)", new File("data/originateCalibrations.nex"));
	final public Input<File> glottoLocationsInput = new Input<>("glottoLocations", "CSV file downloaded from GlottoLog with locations", new File("data/languages_and_dialects_geo.csv"));
	final public Input<File> customLocationsInput = new Input<>("customLocations", "CSV file with locations not in GlottoLog file", new File("data/languages_and_dialects_geo2.csv"));
	final public Input<XMLFile> templateInput = new Input<>("template","template of BEAST XML for EDGE analysis", new XMLFile("data/template.xml"));
	final public Input<OutFile> outputInput = new Input<>("out","output directory, edge.xml + other outputs for sanity checking are written here", new OutFile("output"));
	final public Input<Boolean> verboseInput = new Input<>("verbose", "if true, generates svg files for all families + outputs more on screen", false);

	
	List<Tree> familyTrees;
	List<Set<String>> families;
	List<String> familyNames;
	List<Integer> familySizes;
	List<Set<String>> clades; // rest of world, sahul, americas + specials for nodes internal to families, e.g., arabic, daghestani
	List<String> cladeNames;
	
	List<MRCAPrior> mrcaCalibrations;
	List<MRCAPrior> origateCalibrations;
	
	Set<String> specials;
	Set<String> austronesianTaxa;
	Set<String> afroasiaticTaxa;
	
	class TaxonLocation {
		public TaxonLocation(String string, double lat, double lon) {
			id = string;
			latitude = lat;
			longitude = lon;
		}
		String id;
		double latitude;
		double longitude;
	}
	
	List<TaxonLocation> taxonLocations;

	@Override
	public void initAndValidate() {
	}
	
	public EdgeXMLGenerator() throws IOException {
	}
	
	@Override
	public void run() throws Exception {
		// 1. select taxa
		familyTrees = processTreeFile();
		Set<String> taxa = getAllTaxonSet();

		// 6. filter out taxa without locations
		processLocations(taxa);

		// 2. define clades
		processClades(taxa);
		
		// 3. generate GlottoLog constraint
		String newick = glottologConstraint(taxa);
		
//		System.out.println("Families:");
//		for (int i = 0; i < familyTrees.size(); i++) {
//			if (familySizes.get(i) > 1) {
//				System.out.println(familyNames.get(i));
//			}
//		}
//		System.out.println("\nSingletons:");
//		for (int i = 0; i < familyTrees.size(); i++) {
//			if (familySizes.get(i) == 1) {
//				System.out.println(familyNames.get(i));
//			}
//		}
//		System.out.println("\nIgnored:");
//		for (int i = 0; i < familyTrees.size(); i++) {
//			if (familySizes.get(i) < 1) {
//				System.out.println(familyNames.get(i));
//			}
//		}
		// 4. add time calibrations on families and clade MRCAs
		mrcaCalibrations = processMRCACalibrations(calibrationsInput.get(), taxa, false);
		
		// 5. add time calibrations on originates of MRCAs
		origateCalibrations = processOriginateCalibrations(originateCalibrationsInput.get());
		
		checkAllTaxaHaveAnOriginateCalibration(taxa, origateCalibrations);
		
		checkForDuplicateCalibrationIDs(mrcaCalibrations, origateCalibrations);
		
		listWithinFamilyCalibrations();
		
		// 7. process template (e.g., add prior on tree stuff)
		processTemplate(taxa, newick);
		
		// 8. run monos
		// 9. run splits

//		FileWriter outfile = new FileWriter("/tmp/x");
//        outfile.write(getAllTaxonSetXML(taxa, "EDGE"));
//        outfile.close();

        Log.warning("Done!");
	}
	
	private void listWithinFamilyCalibrations() {
		System.out.println("List of clades encapsulated in other clades");
		for (int i = 0; i < clades.size(); i++) { 
			Set<String> c = clades.get(i);
			if (!(cladeNames.get(i).equals("Americas") || cladeNames.get(i).equals("RestOfWorld"))) {
				for (int j = 0; j < clades.size(); j++) { 
					Set<String> c2 = clades.get(j);
					if (c2.size() < c.size() && c.containsAll(c2)) {
						System.out.println(cladeNames.get(i) + " << " + cladeNames.get(j));
					}
				}			
				for (int j = 0; j < families.size(); j++) { 
					Set<String> c2 = families.get(j);
					if (c2.size() > c.size() && c2.containsAll(c)) {
						System.out.println(familyNames.get(j) + " << " + cladeNames.get(i));
					}
				}		
			}
		}
		System.out.println("End of list of clades encapsulated in other clades");
		
	}

	private void checkForDuplicateCalibrationIDs(List<MRCAPrior> mrcaCalibrations2,
			List<MRCAPrior> origateCalibrations2) {
		// check IDs are unique
		Set<String> ids = new LinkedHashSet<>();
		for (MRCAPrior p : mrcaCalibrations2) {
			if (ids.contains(p.getID())) {
				Log.warning("Duplicate calibration id: " + p.getID());
			}
			ids.add(p.getID());
		}
		for (MRCAPrior p : origateCalibrations2) {
			if (ids.contains(p.getID())) {
				Log.warning("Duplicate calibration id: " + p.getID());
			}
			ids.add(p.getID());
		}
		
	}

	private void checkAllTaxaHaveAnOriginateCalibration(Set<String> taxa, List<MRCAPrior> origateCalibrations) {
		Set<String> taxaCopy = new LinkedHashSet<>();
		taxaCopy.addAll(taxa);
		
		Set<String> taxaWithCalibration = new LinkedHashSet<>();
		for (MRCAPrior p : origateCalibrations) {
			taxaWithCalibration.addAll(p.taxonsetInput.get().asStringList());
		}
		taxaCopy.removeAll(taxaWithCalibration);
		for (String taxon: taxaCopy) {
			Log.warning("WARNING: Taxon " + taxon + " has no originate prior");
		}
	}

	private void processTemplate(Set<String> taxa, String glottologConstraints) throws IOException {
        BufferedReader fin = new BufferedReader(new FileReader(templateInput.get()));
        PrintStream out = new PrintStream(outputInput.get().getPath() + "/edge.xml");
        String allTaxa = getAllTaxonSetXML(taxa, "EDGE");
        
        String taxonSets = createTaxonSetsXML();
        String MRCACalibrations = createCalbrationXML(mrcaCalibrations, false, false);
        String originateCalibrations = createCalbrationXML(this.origateCalibrations, true, false);
        String locations = listLocations();
        String rootPrior = createCalbrationXML(mrcaCalibrations, false, true);
        String excludedForMultiYule = getExcludedForMultiYule();
        
        while (fin.ready()) {
            String str = fin.readLine();
            if (str.indexOf('%') >= 0) {
            	str = replaceAll(str, "%all-taxa%", allTaxa);
            	str = replaceAll(str, "%taxon-sets%", taxonSets);
            	str = replaceAll(str, "%glottologconstraints%", glottologConstraints);
            	str = replaceAll(str, "%MRCACalibrations%",MRCACalibrations); 
            	str = replaceAll(str, "%originateCalibrations%", originateCalibrations); 
            	str = replaceAll(str, "%locations%", locations);
            	str = replaceAll(str, "%root_prior%", rootPrior);
            	str = replaceAll(str, "%excluded_for_multiYule%", excludedForMultiYule);
            	
            }
        	out.println(str);
        }
        out.close();
		fin.close();
	}

	private String getExcludedForMultiYule() {
		StringBuilder b = new StringBuilder();
		for (String str : specials) {
			b.append("                  <exclude idref=\"" + str + ".prior\"/>\n");
		}
		return b.toString();
	}

	private String replaceAll(String str, String pattern, String replacement) {
		int i = str.indexOf(pattern);
		if (i >= 0) {
			str = (i>0?str.substring(0, i) : "") + 
					replacement +
					str.substring(i + pattern.length());
		}
		return str;
	}

	private String createCalbrationXML(List<MRCAPrior> calibrations, boolean useOriginate, boolean rootOnly) {
		StringBuilder b = new StringBuilder();
		for (MRCAPrior p : calibrations) {
			if ((p.taxonsetInput.get().getID().equals("EDGE") && rootOnly) || 
				(!p.taxonsetInput.get().getID().equals("EDGE") && !rootOnly)){
				String id = p.taxonsetInput.get().getID();
		        b.append("<distribution id=\"" + id + (useOriginate ? ".originate" : "") +".prior\" spec=\"beast.math.distributions.MRCAPrior\" monophyletic=\"true\" tree=\"@Tree.t:EDGE\""
		        		+ (useOriginate ? " useOriginate=\"true\" " : "")
		        		+ ">\n"); 
		        b.append("<taxonset idref=\""+id+"\"/>\n");
		        ParametricDistribution distr = p.distInput.get();
		        if (distr != null) {
		        	if (distr instanceof Normal) {
		        		Normal ln = (Normal) distr;
		    	        b.append("<Normal name=\"distr\" "
		    	        		+ "mean=\"" + ln.meanInput.get().getValue() + "\" "
		    	        		+ "sigma=\"" + ln.sigmaInput.get().getValue() + "\"/>\n");
		        	} else if (distr instanceof LogNormalDistributionModel) {
		        		LogNormalDistributionModel ln = (LogNormalDistributionModel) distr;
		    	        b.append("<LogNormal name=\"distr\" meanInRealSpace=\"" + ln.hasMeanInRealSpaceInput.get() + "\" "
		    	        		+ "M=\"" + ln.MParameterInput.get().getValue() + "\" "
		    	        		+ "S=\"" + ln.SParameterInput.get().getValue() + "\"/>\n");
		        		
		        	} else if (distr instanceof Uniform) {
		        		Uniform ln = (Uniform) distr;
		    	        b.append("<Uniform name=\"distr\" "
		    	        		+ "upper=\"" + ln.upperInput.get() + "\" "
		    	        		+ "lower=\"" + ln.lowerInput.get() + "\"/>\n");
		        	} else if (distr instanceof Exponential) {
		        		Exponential ln = (Exponential) distr;
		    	        b.append("<Exponential name=\"distr\" "
		    	        		+ "mean=\"" + ln.lambdaInput.get().getValue() + "\"/>\n");
		        	}
		        }
		        b.append("</distribution>\n");
			}
		}
		return b.toString();
	}

	private String createTaxonSetsXML() {
		StringBuilder b = new StringBuilder();
		for (int i = 0; i < familyNames.size(); i++) {
			if (familySizes.get(i) > 0) {
				b.append("<taxonset id=\"" + familyNames.get(i)+ "\"><plate var=\"n\" range=\"");
				for (String s : families.get(i)) {
					b.append(s);
					b.append(',');
				}
				b.deleteCharAt(b.length() - 1);
				b.append("\"><taxon idref=\"$(n)\"/></plate></taxonset>\n");
			}
		}

		for (int i = 0; i < clades.size(); i++) {
			b.append("<taxonset id=\"" + cladeNames.get(i)+ "\"><plate var=\"n\" range=\"");
			for (String s : clades.get(i)) {
				b.append(s);
				b.append(',');
			}
			b.deleteCharAt(b.length() - 1);
			b.append("\"><taxon idref=\"$(n)\"/></plate></taxonset>\n");
		}		
		return b.toString();
	}

	private String listLocations() {
		StringBuilder b = new StringBuilder();
		for (TaxonLocation t : taxonLocations) {
			b.append(t.id + " = " + t.latitude + " " + t.longitude + ",\n");
		}		
		return b.toString();
	}

	private void processLocations(Set<String> taxa) throws IOException {
		taxonLocations = new ArrayList<>();
		Set<String> taxa2 = new LinkedHashSet<>();
		taxa2.addAll(taxa);
		
		
		getLocations(glottoLocationsInput.get(), taxa2);
		getLocations(customLocationsInput.get(), taxa2);

		// sanity check: all locations assigned?
		PrintStream out = new PrintStream(outputInput.get().getPath() + "/removedDueToLackOfLocation.log");
		int k = 0;
		for (String t : taxa2) {
			Log.warning("Removed due to lack of location: " + t);
			taxa.remove(t);
			out.println(t);
			k++;
		}
		if (k > 0) {
			Log.warning(k + " languages removed due to lack of location");
			out.println(k + " languages removed due to lack of location");
		}
		out.close();
	}

	private void getLocations(File file, Set<String> taxa) throws IOException {
        BufferedReader fin = new BufferedReader(new FileReader(file));
        // eat up header line
        String str = fin.readLine();
        while (fin.ready()) {
            str = fin.readLine();
            if (!str.startsWith("#")) {
            	String [] strs = str.split(",");
            	if (taxa.contains(strs[0]) && str.lastIndexOf(',') != str.length()-1) {
            		double lat = Double.parseDouble(strs[strs.length - 2]);
            		double lon = Double.parseDouble(strs[strs.length - 1]);
            		TaxonLocation t = new TaxonLocation(strs[0], lat, lon);
            		taxonLocations.add(t);
            		taxa.remove(strs[0]);
            		// System.out.print(t.id + ",");
            	}            	
            }
        }
        fin.close();
	}

	private List<MRCAPrior> processOriginateCalibrations(File nexusFile) throws IOException {
		List<MRCAPrior> originateCalibrations = processMRCACalibrations(nexusFile, null, true);
		return originateCalibrations;
	}

	private List<MRCAPrior> processMRCACalibrations(File nexusFile, Set<String> alltaxa, boolean useOriginate) throws IOException {
		List<MRCAPrior> calibrations = new ArrayList<>();
		NexusParser parser = new NexusParser();
		
        BufferedReader fin = new BufferedReader(new FileReader(nexusFile));
        while (fin.ready()) {
            String str = fin.readLine();

            if (!str.startsWith("#") && str.toLowerCase().matches("^\\s*calibrate\\s.*")) {
		    	// define calibration represented by an MRCAPRior, 
		    	// taxon sets need to be specified earlier, but can also be a single taxon
		    	// e.g.
		    	// begin assumptions;
		    	// calibrate germanic = normal(1000,50)
		    	// calibrate hittite = normal(3450,100)
		    	// calibrate english = fixed(0)
		    	// end;
		    	String [] strs = str.split("=");
		    	if (strs.length > 1) {
		    		String str0 = strs[0].trim();
		    		String [] strs2 = str0.split("\\s+");
		    		if (strs2.length != 2) {
		    			throw new RuntimeException("expected 'calibrate <name> = ...' but did not get two words before the = sign: " + str);
		    		}
		    		// first, get the taxon
		    		String taxonSetName = strs2[1].replaceAll("'\"", "");
		    		int i = indexof(familyNames, taxonSetName);
		    		Set<String> taxa = null;
		    		if (i < 0) {
		    			i = indexof(cladeNames, taxonSetName);
		    			if (i < 0) {
		    				if (taxonSetName.equals("EDGE")) {
		    					taxa = alltaxa;
		    				} else {
		    					Log.warning("Could not find clade " + taxonSetName + " for calibration");
		    				}
		    			} else {
		    				taxa = clades.get(i);
		    			}
		    		} else {
		    			taxa = families.get(i); 
		    		}
		    		
		    		if (taxa != null && taxa.size() > 0) {
		    			TaxonSet taxonset = createTaxonSet(taxa, taxonSetName);
			    		if (taxonset == null) {
			    			throw new RuntimeException("Could not find taxon/taxonset " + taxonSetName + " in calibration: " + str);
			    		}
			    		
			    		// next get the calibration
			    		str0 = strs[strs.length - 1].trim();
			    		String [] strs3 = str0.split("[\\(,\\)]");
			
			    		try {
			                MRCAPrior prior = parser.getMRCAPrior(taxonset, strs3, useOriginate);
			
			                // should set Tree before initialising, but we do not know the tree yet...
			                if (calibrations == null) {
			                    calibrations = new ArrayList<>();
			                }
			                calibrations.add(prior);
			            } catch (RuntimeException ex) {
			                throw new RuntimeException(ex.getMessage() + "in calibration: " + str);
			            }
		    		}
		    	}
	        }
        }
	    fin.close();

        return calibrations;
	}

	private String glottologConstraint(Set<String> taxa) throws IOException {
		StringBuilder b = new StringBuilder();
		b.append("(");
		Set<String> taxa2 = new LinkedHashSet<>();
		taxa2.addAll(taxa);
		familySizes = new ArrayList<>();
		PrintStream out = new PrintStream(outputInput.get() + "/excludedFamilies.log");
		out.println("The list below should be equal to that in " + familyExclusionsInput.get());
		for (Tree tree : familyTrees) {
			int [] count = new int [1];
			if (!toNewick(tree.getRoot(), taxa2, b, count)) {
				String name = tree.getRoot().getID();
				String glottocode = name.substring(name.indexOf('[') + 1, name.indexOf(']'));
				name = name.substring(1, name.indexOf('['));
				out.println(glottocode + "," + name);
			}
			b.append(',');
			familySizes.add(count[0]);
		}		
		b.deleteCharAt(b.length() - 1);
		b.append(")");
		out.close();

		// sanity check
		for (String s : taxa2) {
			Log.warning("Taxon " + s + " is not in glottologConstraint");
		}
		
		// clean up commas 
		String newick = b.toString();
		while (newick.contains(",,")) {
			newick = newick.replaceAll(",,", ",");
		}
		return newick;
	}
	
	/** convert a family tree to newick for glottolog constraints, taking only
	 * taxa in account. Keeps track of how many leafs end up in the newick tree
	 * through the count variable.
	 */
	private boolean toNewick(Node node, Set<String> taxa, StringBuilder b, int [] count) {
		String id = node.getID();
		if (id != null) {
			id = id.substring(id.indexOf('[') + 1, id.indexOf(']'));
			if (taxa.contains(id)) {
				taxa.remove(id);
				b.append(id);
				count[0]++;
				return true;
			}
		}
		StringBuilder [] cb = new StringBuilder[node.getChildCount()];
		
		int k = 0;
		for (int i = 0; i < node.getChildCount(); i++) {
			cb[i] = new StringBuilder();
			if (toNewick(node.getChild(i), taxa, cb[i], count)) {
				k++;
			} else {
				cb[i] = null;
			}
		}
		
		
		if (k > 0) {
			if (k > 1) {
				b.append('(');
				for (StringBuilder sb : cb) {
					if (sb != null) {
						b.append(sb);
						b.append(',');
					}
				}
				b.deleteCharAt(b.length() - 1);
				b.append(')');
			} else {
				for (StringBuilder sb : cb) {
					if (sb != null) {
						b.append(sb);
					}
				}
			}
			return true;
		}
		return false;
	}

	private void processClades(Set<String> alltaxa)  throws IOException {
		createFamilies(alltaxa);
		loadClades(alltaxa);
		loadSpecials(alltaxa);
	}

	private void loadSpecials(Set<String> alltaxa) throws IOException {
		specials = new LinkedHashSet<>();
		File specialsFile = new File(cladesInput.get().getPath() + "/specials.csv");
		for (String str : BeautiDoc.load(specialsFile).split("\n")) {
			if (!str.startsWith("#")) {
				String [] strs = str.split(",");
				String name = strs[0];
				String glottocode = strs[1];
				Node node = find(glottocode);
				Set<String> taxa = new LinkedHashSet<>();
				getTaxa(node, taxa);
				taxa.retainAll(alltaxa);
				cladeNames.add(name);
				clades.add(taxa);
				this.specials.add(name);
			}
		}		
	}

	private Node find(String glottocode) {
		for (Tree tree : familyTrees) {
			Node result = find(tree.getRoot(), glottocode);
			if (result != null) {
				return result;
			}
		}
		return null;
	}

	private Node find(Node node, String glottocode) {
		String id = node.getID();
		if (id != null && id.contains(glottocode)) {
			return node;
		}
		for (Node child : node.getChildren()) {
			Node found = find(child, glottocode);
			if (found != null) {
				return found;
			}
		}
		return null;
	}

	private void loadClades(Set<String> alltaxa) throws IOException {
		clades = new ArrayList<>();
		cladeNames = new ArrayList<>();
		for (String f : cladesInput.get().list()) {
			if (f.endsWith("fam")) {
				String name = f;
				if (name.contains(".")) {
					name = name.substring(0, name.lastIndexOf("."));
				}
				if (name.contains("\\")) {
					name = name.substring(name.lastIndexOf("//"));
				}
				if (name.contains("/")) {
					name = name.substring(name.lastIndexOf("/"));
				}
				f = cladesInput.get().getPath() + "/" + f;
				cladeNames.add(name);
				
				Set<String> taxa = new LinkedHashSet<>();
				boolean allfound = true;
				for (String family : BeautiDoc.load(new File(f)).split("\n")) {
					int i = indexof(familyNames, family);
					if (i < 0) {
						if (allfound) {
							Log.warning("Trying to define clade " + name + " from families in file " + f);
							allfound = false;
						}
						Log.warning("Cannot find family " + family + " for clade " + name);
					} else {
						taxa.addAll(families.get(i));
					}					
				}
				taxa.retainAll(alltaxa);
				clades.add(taxa);
				dotsOnMap(taxa, name);
				
				// check dots without Austronesia
				Set<String> noaustaxa = new LinkedHashSet<>();
				noaustaxa.addAll(taxa);
				noaustaxa.removeAll(austronesianTaxa);
				dotsOnMap(noaustaxa, name+"-withoutAustronesia");
				noaustaxa.removeAll(afroasiaticTaxa);
				dotsOnMap(noaustaxa, name+"-withoutAfroAsiatic");
				
				if (name.equals("RestOfWorld")) {
					Set<String> noROW = new LinkedHashSet<>();
					noROW.addAll(alltaxa);
					noROW.removeAll(taxa);
					dotsOnMap(noROW, "Africa");
				}
			}			
		}
	}

	private void dotsOnMap(Set<String> taxa, String name) throws IOException {
		PrintStream out = new PrintStream(outputInput.get().getPath() + "/" + name + ".svg");
		int w = 1000, h = 500;
		out.println("<svg width=\""+w+"\" height=\""+h+"\"\n" + 
				  "xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n"+       
				  "<image xlink:href=\"https://www.cs.auckland.ac.nz/~remco/geo/World98.png\" height=\""+h+"\" width=\""+w+"\"/>\n");
		for (String t : taxa) {
			TaxonLocation loc = null;
			for (int i = 0; i < taxonLocations.size(); i++) {
				if (taxonLocations.get(i).id.equals(t)) {
					loc = taxonLocations.get(i);
					break;
				}
			}
			double x = (loc.longitude + 180.0)*w/360.0;
			double y = (90.0-loc.latitude)*h/180.0;
			out.println("<circle cx=\"" + x + "\" cy=\"" + y + "\" id=\""+t+"\" r=\"1.5\" style=\"fill:#00f\"/>");
		}
		out.println("</svg>\n");
		out.close();
	}

	private int indexof(List<String> list, String str) {
		for (int i = 0; i < list.size(); i++) {
			if (list.get(i).equals(str)) {
				return i;
			}
		}
		return -1;
	}

	private void createFamilies(Set<String> alltaxa) throws IOException {
		families = new ArrayList<>();
		familyNames = new ArrayList<>();
		for (Tree tree : familyTrees) {
			String name = tree.getRoot().getID();
			name = name.substring(1, name.indexOf('[') - 1);
			name = name.replaceAll(" ", "_");
			// name = Normalizer.normalize(name, Normalizer.Form.NFD);
			name = name.replaceAll("[áã]", "a");
			name = name.replaceAll("í", "i");
			name = name.replaceAll("[ôó]", "o");
			name = name.replaceAll("[éê]", "e");
			name = name.replaceAll("[úü]", "u");
			
			familyNames.add(name);
			Set<String> taxa = new LinkedHashSet<>();
			getTaxa(tree.getRoot(), taxa);
			taxa.retainAll(alltaxa);
			families.add(taxa);
			if (name.equals("Austronesian")) {
				austronesianTaxa = taxa;
			}
			if (name.equals("Afro-Asiatic")) {
				afroasiaticTaxa = taxa;
			}
			
			if (verboseInput.get()) {
				dotsOnMap(taxa, name);
			}
		}		
	}

	/** processes glottolog newick tree file 
	 * 
	 * **/	
	private List<Tree> processTreeFile() throws IOException {
		List<String> exclusions = new ArrayList<>();
		for (String str : BeautiDoc.load(familyExclusionsInput.get()).split("\n")) {
			if (!str.startsWith("#")) {
				str = str.split(",")[0];
				exclusions.add(str);
			}
		}
		
		List<Tree> familyTrees = new ArrayList<>();
        BufferedReader fin = new BufferedReader(new FileReader(treeFileInput.get()));
        while (fin.ready()) {
            String newick = fin.readLine();
            processFamily(newick, familyTrees, exclusions);
        }
        fin.close();
        return familyTrees;
	}

	/** processes one language family from glottolog newick tree file **/
	private void processFamily(String newick, List<Tree> familyTrees, List<String> exclusions) {
		newick = newick.replaceAll("''", " " );
		TreeParser parser = new TreeParser(newick, false, true, true, 0, false);
		
		String familyName = parser.getRoot().getID();
		for (String glottocode : exclusions) {
			if (familyName.contains(glottocode)) {
				return;
			}
		}
		familyTrees.add(parser);
	}

	public String getAllTaxonSetXML(Set<String> taxa, String id) throws IOException {
		StringBuilder b = new StringBuilder();
		b.append("<taxonset id=\"" + id + "\">\n");
		for (String tid : taxa) {
			b.append("<taxon id=\"" + tid + "\"/>\n");
		}
		b.append("</taxonset>");
		return b.toString();
	}
	
	public Set<String> getAllTaxonSet() throws IOException {
		// first, collect all languages from glottolog tree
		Set<String> taxa = new LinkedHashSet<>();
		for (Tree tree : familyTrees) {
			getTaxa(tree.getRoot(), taxa);
		}
		
		// remove ancient languages, colonies
		String [] exclusions = BeautiDoc.load(languageExclusionsInput.get()).split("\n");
		int removed = 0, notfound = 0;
		boolean allfound = true;
		PrintStream out = new PrintStream(outputInput.get().getPath() + "/languageExclusion.log");
		String reason = "";
		for (String glottocode : exclusions) {
			if (!glottocode.startsWith("#")) {
				if (taxa.contains(glottocode)) {
					taxa.remove(glottocode);
					Log.warning("Removing " + glottocode);
					out.println(glottocode + " "  + reason);
					removed++;
				} else {
					if (allfound) {
						Log.warning("Processing taxa from exclusion list to remove (" + languageExclusionsInput.get().getPath() + ")");
						allfound = false;
					}
					Log.warning("Could not find glottocode " + glottocode + " to remove from langauges");
					out.println(glottocode);
					notfound++;
				}
			} else {
				reason = glottocode; 
			}
		}
		Log.warning("Exclusions: " + removed + " removed " + (notfound > 0 ? notfound + " not found" : ""));
		out.println("Exclusions: " + removed + " removed " + (notfound > 0 ? notfound + " not found (listed above)" : ""));
		out.println();
		out.close();
		return taxa;		
	}
	
	private void getTaxa(Node node, Set<String> taxa) {
		if (node.getID() != null && node.getID().contains("-l-")) {
			String id = node.getID();
			id = id.substring(id.indexOf('[')+1 , id.indexOf(']'));
			taxa.add(id);
		} else {
			for (Node child : node.getChildren()) {
				getTaxa(child, taxa);
			}
		}
	}

	// to keep track of taxa, so no duplicates are created
	List<Taxon> taxa = new ArrayList<>();
	
	private TaxonSet createTaxonSet(Set<String> codes, String id) {
		List<Taxon> taxaAsList = new ArrayList<>();
		for (String s : codes) {
			boolean found = false;
			for (Taxon t : taxa) {
				if (t.getID().equals(t)) {
					taxaAsList.add(t);
					found = true;
					break;
				}
			}
			if (!found) {
				Taxon t = new Taxon(s);
				taxa.add(t);
				taxaAsList.add(t);
			}
		}
		TaxonSet taxonset = new TaxonSet(taxaAsList);
		taxonset.setID(id);
		return taxonset;
	}
	
	public static void main(String[] args) throws Exception {
		new Application(new EdgeXMLGenerator(), "Generate EDGE XML", args);
	}

}
