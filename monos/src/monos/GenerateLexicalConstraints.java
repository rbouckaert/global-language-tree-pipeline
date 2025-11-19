package monos;



import java.io.*;
import java.util.*;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import beastfx.app.util.XMLFile;
import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.core.Log;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.CladeSet;
import beast.base.evolution.tree.Node;
import beastlabs.evolution.tree.PrunedTree;
import beast.base.evolution.tree.Tree;
import almostbeast.math.distributions.AlmostMRCAPrior;
import almostbeast.math.distributions.AlmostMultiMRCAPriors;
import beast.base.evolution.tree.MRCAPrior;
import beast.base.inference.distribution.Normal;
import beast.base.evolution.tree.TreeParser;
import beast.base.parser.XMLParser;
import beast.base.parser.XMLParserException;
import beast.base.parser.XMLProducer;

@Description("java implementation of the monos package -- add lexical constraints to XML based on "
		+ "lexical analysis (if compatible with glottolog constraints)")
public class GenerateLexicalConstraints extends beast.base.inference.Runnable {
	final public Input<TreeFile> glottologTreeInput = new Input<>("glottologTree", 
			"Newick file containing the monophyletic constraints from Glottolog", Validate.REQUIRED);


	/** input required for finding out the taxon set **/ 
	final public Input<XMLFile> modelInput = new Input<>("xml",
			"file name of BEAST XML file containing the model for which to create constraints for",
			new XMLFile("Glottolog.xml"), Validate.REQUIRED);
	final public Input<String> treeIDInput = new Input<>("treeID", "name of the tree in the XML file that contains all taxa. "
			+ "This taxon set should contain all taxa in the target output.", Validate.REQUIRED);

	/** inputs required for finding lexical analyses **/ 
	final public Input<File> treeConfigInput = new Input<>("treeConfig", 
			"configuration file that identifies lexical analyses. "
			+ "This should be a tab delimited file with three columns: name of the posterior tree set, "
			+ "name of the taxon map file, and scale factor for trees.", Validate.REQUIRED);
	// used to suppress warning messages for tree set cfg file checking only
	final public Input<File> languageExclusionsInput = new Input<>("languageExclusions", "text file with GlottoLog codes of languages that should not be included in the analysis");

	
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees "
			+ "to used as burn-in (and will be ignored)", 10);

	final public Input<OutFile> outputInput = new Input<>("out", "output file. Print to stdout if not specified");
	
	final public Input<Double> thresHoldInput = new Input<>("threshold", "threshold above which clade support in lexical analyses "
			+ "is deemed high enough to justify adding a monophyly constraint in the pruned tree.", 0.95);
	
	MCMC mcmc;
	Tree tree;
	PrunedTree prunedTree;
	
	
	int treeCount;
	Map<String, Taxon> taxonMap = new LinkedHashMap<>();
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		// make sure all files exist
		checkArguments();

		Tree glottologTree = getGlottologTree();
		checkIsosInLexicalAnalyses(glottologTree);

		/** set of taxa in the XML **/
		TaxonSet allTaxaInAnalysis = getTaxonset();
		
		/** set of taxa in lexical analyses that are also in allTaxaInAnalysis **/
		List<Taxon> prunedTaxa = getPosteriorTaxonset(allTaxaInAnalysis);		
		createPrunedTree(prunedTaxa);
		
		processPosteriors(glottologTree);
		
		output();
		Log.warning("Done");
		Log.warning("\nYou need to manually edit KML constraints back in if you have KMLRegions");
		
	}


	private void checkIsosInLexicalAnalyses(Tree glottologTree) throws IOException {		
		String cfg = BeautiDoc.load(treeConfigInput.get());
		String [] strs = cfg.split("\n");
		Log.warning("Checking tree config file codes in " + treeConfigInput.get().getPath());
		
		
		Set<String> tabu = new LinkedHashSet<>();
		tabu.add("glottocode");
		if (languageExclusionsInput.get() != null) {
			String str = BeautiDoc.load(languageExclusionsInput.get());
			String [] strs2 = str.split("\n");
			for (String s : strs2) {
				if (!s.startsWith("#") && s.trim().length() > 0) {
					tabu.add(s.trim());
				}
			}
		}
		
		
		List<String> taxa = glottologTree.getTaxonset().asStringList();
		List<String> missing = new ArrayList<>();
		for (String str : strs) {
			if (str.trim().length() > 0 && !str.startsWith("#")) {
				String isomapfile = str.split("\t")[1];
				str = BeautiDoc.load(new File(isomapfile));
				for (String s : str.split("\n")) {
					String iso = s.split("\t")[0];
					if (taxa.indexOf(iso) < 0 && !tabu.contains(iso)) {
						Log.warning("code " + iso + " from " + isomapfile + " is not in glottolog tree");
						missing.add(iso);						
					}
				}
			}
		}
		if (missing.size() > 0) {
			Log.warning("Missing codes do not necessarily mean something is wrong (e.g., can include ancient languages). "
					+ "This list is just useful for checking: " + missing);
		}
		Log.warning("Config file check done");
	}

	private void checkArguments() throws IOException {
		fileCheck(glottologTreeInput.get());
		fileCheck(modelInput.get());
		
		fileCheck(treeConfigInput.get());
		
		String cfg = BeautiDoc.load(treeConfigInput.get());
		String [] strs = cfg.split("\n");
		Log.warning("Checking tree config file " + treeConfigInput.get().getPath());
		for (String str : strs) {
			if (str.trim().length() > 0 && !str.trim().startsWith("#")) {
				String treefile = str.split("\t")[0];
				fileCheck(new File(treefile));
				String isomapfile = str.split("\t")[1];
				fileCheck(new File(isomapfile));
			}
		}
	}

	private void fileCheck(File file) {
		if (!file.exists()) {
			throw new IllegalArgumentException("File >>" + file.getPath() + "<< does not exist");
		}
		if (!file.canRead()) {
			throw new IllegalArgumentException("File >>" + file.getPath() + "<< exists but cannot be read");
		}
	}

	/** extract clade sets from posterior tree files
	 *  clades with >95% support -- if compatible with glottologTree -- are added as constraints
	 */
	private void processPosteriors(Tree glottologTree) throws IOException {
		
		BEASTInterface o2 = getObjectWithID(mcmc, "prior");;
		if (o2 == null) {
			throw new IllegalArgumentException("could not find tree with id=\"prior\"");
		}
		if (!(o2 instanceof CompoundDistribution)) {
			throw new IllegalArgumentException("Object identified by id=\"prior\" is not a CompoundDistribution");
		}
		CompoundDistribution prior = (CompoundDistribution) o2;
		AlmostMultiMRCAPriors extraPriors = new AlmostMultiMRCAPriors();
		extraPriors.setID("lexical-priors");
		extraPriors.treeInput.setValue(prunedTree, extraPriors);
		// if 'newick' were not a required input in MultiMRCAPrior
		// we could skip next line, which just adds a dummy constraint
		extraPriors.newickInput.setValue("()", extraPriors);
		prior.pDistributions.get().add(extraPriors);
		Logger tracelog = (Logger) getObjectWithID(mcmc, "tracelog");
		if (tracelog != null) {
			tracelog.loggersInput.get().add(extraPriors);
		}
		
		String cfg = BeautiDoc.load(treeConfigInput.get());
		String [] strs = cfg.split("\n");
		int firstTaxonSetCount = 0;
		for (String str : strs) {
			if (str.trim().length() > 0 && !str.trim().startsWith("#")) {
				String treefile = str.split("\t")[0];
				String isomapfile = str.split("\t")[1];
				String scaleFactorStr = str.split("\t")[2];
				double scaleFactor = Double.parseDouble(scaleFactorStr);
				
				
				List<MRCAPrior> ageConstraints = new ArrayList<>();
				Map<String,String> taxonMap = getIsoMap(isomapfile);
				CladeSet c = getCladeSet(treefile, taxonMap);
	
				// determine monophyly constraints on Pruned tree
				double thresHold = thresHoldInput.get();
				for (int i = 0; i < c.getCladeCount(); i++) {
					int freq = c.getFrequency(i);
					if (freq >= thresHold * treeCount) {
						String clade = c.getClade(i);
						TaxonSet taxonset = cladeToTaxonSet(clade);
						if (taxonset.getTaxonCount() > 1) {
							if (isCompatible(taxonset, glottologTree)) {
								MRCAPrior mrcaPrior = new AlmostMRCAPrior();
								mrcaPrior.initByName("tree", prunedTree, 
										"taxonset", taxonset, 
										"monophyletic", true);
								if (treefile.indexOf('/') > 0) {
									mrcaPrior.setID("lexical_" + treefile.substring(0,treefile.indexOf('/')) + "_");
								} else {
									mrcaPrior.setID("lexical_prior_");
								}
								
								if (firstTaxonSetCount < taxonset.getTaxonCount()) {
									// make sure biggest taxonset goes first
									// this should help with monos/filter/filterg.pl to generate valid XML
									firstTaxonSetCount = taxonset.getTaxonCount();
									extraPriors.calibrationsInput.get().add(0, mrcaPrior);
								} else {
									extraPriors.calibrationsInput.get().add(mrcaPrior);
								}
								//prior.pDistributions.get().add(mrcaPrior);
								ageConstraints.add(mrcaPrior);
							} else {
//								Log.debug("Incompatible: " + freq + " " + taxonset.asStringList());
							}
						}					
					}
				}
				
				// determine age distributions
				getDistributions(treefile, ageConstraints, taxonMap, scaleFactor);
			}
		}		
	}


	private Map<String, String> getIsoMap(String isomapfile) throws IOException {
		String str = BeautiDoc.load(isomapfile);
		String [] strs = str.split("\n");
		Map<String, String> isoMap = new LinkedHashMap<>();
		for (String str2 : strs) {
			String [] strs2 = str2.split("\t");
			if (strs2.length > 1) {
				isoMap.put(strs2[1], strs2[0]);
			}
		}
		return isoMap;
	}

	private TaxonSet cladeToTaxonSet(String clade) {
		TaxonSet taxonset = new TaxonSet();
		String [] strs = clade.split(",");
		for (String str : strs) {
			str = str.trim();
			if (str.startsWith("{")) {
				str = str.substring(1);
			}
			if (str.endsWith("}")) {
				str = str.substring(0, str.length() - 1);
			}
			taxonset.taxonsetInput.get().add(taxonMap.get(str));
		}		
		taxonset.initAndValidate();
		return taxonset;
	}

	/** clade with taxonset leafs of prunedTree taxa is compatible with
	 *  glottologTree if the mrca of taxonset in glottologTree contains 
	 *  no other taxa in prunedTree
	 */
	private boolean isCompatible(TaxonSet taxonset, Tree glottologTree) {
		Node mrca = getCommonAncestor(taxonset, glottologTree);
		
		String [] prunedTaxa = prunedTree.getTaxaNames();
		for (Node child : mrca.getChildren()) {
			if (nodesTraversed[child.getNr()]) {
				for (Node leaf : child.getAllLeafNodes()) {
					String id = leaf.getID();
					if (contains(prunedTaxa, id) && !contains(taxonset, id)) {
						// System.err.println("incompatible id (" + id+ "): " + taxonset);
						return false;
					}
				}
			}
		}		

		return true;
	}
	
    private boolean contains(TaxonSet taxonset, String id) {
		for (Taxon t : taxonset.taxonsetInput.get()) {
			if (t.getID().equals(id)) {
				return true;
			}
		}
		return false;
	}

	private boolean contains(String[] prunedTaxa, String id) {
		for (String p : prunedTaxa) {
			if (p.equals(id)) {
				return true;
			}
		}
		return false;
	}

	/** keep track which nodes are already traversed when finding MRCA of set of nodes **/ 
	boolean [] nodesTraversed;

    private Node getCommonAncestor(TaxonSet taxonSet, Tree tree) {
        nodesTraversed = new boolean[tree.getNodeCount()];
    	List<String> taxa = taxonSet.asStringList();
        final List<String> taxaNames = new ArrayList<>();
        for (final String taxon : tree.getTaxaNames()) {
            taxaNames.add(taxon);
        }
        int i = taxaNames.indexOf(taxa.get(0));
        if (i < 0) {
        	Log.warning("Encountered unknown taxon: " + taxa.get(0));
        	Log.warning("Make sure the taxon exists in the glottolog tree");
        }
        Node cur = tree.getNode(i);
        for (int k = 1; k < taxa.size(); ++k) {
        	i = taxaNames.indexOf(taxa.get(k));
            if (i < 0) {
            	Log.warning("Encountered unknown taxon: " + taxa.get(k));
            	Log.warning("Make sure the taxon exists in the glottolog tree");
            }
            cur = getCommonAncestor(cur, tree.getNode(i));
        }
        return cur;
    }

    private Node getCommonAncestor(TaxonSet taxonSet, Tree tree, Map<String,String> taxonMap) {
        nodesTraversed = new boolean[tree.getNodeCount()];
    	List<String> taxa = taxonSet.asStringList();
        final List<String> taxaNames = new ArrayList<>();
        for (final String taxon : tree.getTaxaNames()) {
            taxaNames.add(taxonMap.get(taxon));
        }

        Node cur = tree.getNode(taxaNames.indexOf(taxa.get(0)));
        for (int k = 1; k < taxa.size(); ++k) {
        	int i = taxaNames.indexOf(taxa.get(k));
            cur = getCommonAncestor(cur, tree.getNode(i));
        }
        return cur;
    }

    
    protected Node getCommonAncestor(Node n1, Node n2) {
        // assert n1.getTree() == n2.getTree();
        if( ! nodesTraversed[n1.getNr()] ) {
            nodesTraversed[n1.getNr()] = true;
        }
        if( ! nodesTraversed[n2.getNr()] ) {
            nodesTraversed[n2.getNr()] = true;
        }
        while (n1 != n2) {
	        double h1 = n1.getHeight();
	        double h2 = n2.getHeight();
	        if ( h1 < h2 ) {
	            n1 = n1.getParent();
	            if( ! nodesTraversed[n1.getNr()] ) {
	                nodesTraversed[n1.getNr()] = true;
	            }
	        } else if( h2 < h1 ) {
	            n2 = n2.getParent();
	            if( ! nodesTraversed[n2.getNr()] ) {
	                nodesTraversed[n2.getNr()] = true;
	            }
	        } else {
	            //zero length branches hell
	            Node n;
	            double b1 = n1.getLength();
	            double b2 = n2.getLength();
	            if( b1 > 0 ) {
	                n = n2;
	            } else { // b1 == 0
	                if( b2 > 0 ) {
	                    n = n1;
	                } else {
	                    // both 0
	                    n = n1;
	                    while( n != null && n != n2 ) {
	                        n = n.getParent();
	                    }
	                    if( n == n2 ) {
	                        // n2 is an ancestor of n1
	                        n = n1;
	                    } else {
	                        // always safe to advance n2
	                        n = n2;
	                    }
	                }
	            }
	            if( n == n1 ) {
                    n = n1 = n.getParent();
                } else {
                    n = n2 = n.getParent();
                }
	            if( ! nodesTraversed[n.getNr()] ) {
	                nodesTraversed[n.getNr()] = true;
	            } 
	        }
        }
        return n1;
    }




	/** produce output to file or stdout (if outputInput is not specified) **/
	private void output() throws FileNotFoundException {
		PrintStream out = outputInput.get() == null ? System.out : new PrintStream(outputInput.get());
		XMLProducer producer = new XMLProducer();
		String xml = producer.toXML(mcmc);
		out.println(xml);
		out.close();
		if (outputInput.get() == null) {
			Log.warning("Output in : " + outputInput.get().getPath());
		}
	}

	/** construct pruned tree from this.tree and prundeTaxa **/
	private void createPrunedTree(List<Taxon> prunedTaxa) {
		TaxonSet prunedTaxonSet = new TaxonSet();
		for (Taxon t : prunedTaxa) {
			prunedTaxonSet.taxonsetInput.get().add(t);
		}
		prunedTaxonSet.initAndValidate();
		prunedTree = new PrunedTree();
		prunedTree.initByName("tree", tree, "taxonset", prunedTaxonSet);
		prunedTree.setID("PrunedTree");
	}

	/** construct taxon set all taxa in posterior tree set **/	
	private List<Taxon> getPosteriorTaxonset(TaxonSet allTaxaInAnalysis) throws IOException {
		String cfg = BeautiDoc.load(treeConfigInput.get());
		String [] strs = cfg.split("\n");
		Set<String> taxa = new LinkedHashSet<>();
		for (String str : strs) {
			if (str.trim().length() > 1 && !str.startsWith("#")) {
				String mapFile = str.split("\t")[1];			
				String map = BeautiDoc.load(mapFile);
				String [] mapStrs = map.split("\n");
				for (String m : mapStrs) {
					taxa.add(m.split("\t")[0]);
				}
			}
		}
		List<Taxon> taxonset = new ArrayList<>();
		for (String t : taxa) {
			for (Taxon taxon : allTaxaInAnalysis.taxonsetInput.get()) {
				if (taxon.getID().equals(t)) {
					taxonset.add(taxon);
				}
			}
		}
		return taxonset;
	}

	/** grab taxon set from the XML file and initialise taxonMap **/
	private TaxonSet getTaxonset() throws SAXException, IOException, ParserConfigurationException, XMLParserException {
		XMLParser parser = new XMLParser();
		mcmc = (MCMC) parser.parseFile(modelInput.get());
		
		BEASTInterface o2 = getObjectWithID(mcmc, treeIDInput.get());;
		if (o2 == null) {
			throw new IllegalArgumentException("could not find tree with id " + treeIDInput.get());
		}
		if (!(o2 instanceof Tree)) {
			throw new IllegalArgumentException("Object identified by id=\"" + treeIDInput.get() +"\" is not a Tree");
		}
		tree = (Tree) o2;
		TaxonSet taxonset = tree.getTaxonset();
		
		for (Taxon t : taxonset.taxonsetInput.get()) {
			taxonMap.put(t.getID(), t);
		}
		return taxonset;
	}

	/** traverse BEAST model in search for object with given id **/
	private BEASTInterface getObjectWithID(BEASTInterface o, String id) {
		for (BEASTInterface o2 : o.listActiveBEASTObjects()) {
			if (o2.getID() != null && o2.getID().equals(id)) {
				return o2;
			} else {
				BEASTInterface t = getObjectWithID(o2, id);
				if (t != null) {
					return t;
				}
			}
		}
		return null;
	}

	/** return tree containing monophyletic constraints based on glottolog **/
	private Tree getGlottologTree() throws IOException {
		String newick = BeautiDoc.load(glottologTreeInput.get().getAbsolutePath());
		TreeParser parser = new TreeParser(newick, true, true, true, 0, false);
		return parser;
	}

	
	/** Approximate age distribution of MRCAs by a normal distribution based on mean & stdDev of
	 *  observed common ancestor heights **/
	private void getDistributions(String path, List<MRCAPrior> ageConstraints, Map<String,String> taxonMap, double scaleFactor) throws IOException {
		List<List<Double>> ages = new ArrayList<>();
		for (int i = 0; i < ageConstraints.size(); i++) {
			ages.add(new ArrayList<>());
		}
		
		Log.warning("Processing " + path);
		
		MemoryFriendlyTreeSet srcTreeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(path, burnInPercentageInput.get());
		srcTreeSet.reset();
		Tree tree = srcTreeSet.next();
		while (srcTreeSet.hasNext()) {
			tree = srcTreeSet.next();
			if (tree == null) {
				break;
			}
			for (int i = 0; i < ageConstraints.size(); i++) {
				TaxonSet taxonset = ageConstraints.get(i).taxonsetInput.get();
				Node mrca = getCommonAncestor(taxonset, tree, taxonMap);
				double age = mrca.getHeight();
				ages.get(i).add(age);
			}
			treeCount++;
		}
		for (int i = 0; i < ageConstraints.size(); i++) {
            List<Double> trace = ages.get(i);
            double sum = 0, sum2 = 0;
            for (double f : trace) {
                sum += f;
                sum2 += f * f;
            }
            double mean = sum / trace.size();
            double stdDev = Math.sqrt(sum2 / trace.size() - mean * mean);
            mean *= scaleFactor;
            stdDev *= scaleFactor;
            Normal distr = new Normal();
            distr.initByName("mean", mean + "", "sigma", stdDev + "");
            ageConstraints.get(i).distInput.setValue(distr, ageConstraints.get(i));
		}		
	}

	
	/** read file from path and decompose into Cladeset **/
	private CladeSet getCladeSet(String path, Map<String,String> taxonMap) throws IOException {
		Log.warning("Processing " + path);
		MemoryFriendlyTreeSet srcTreeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(path, burnInPercentageInput.get());
		srcTreeSet.reset();
		Tree tree = srcTreeSet.next();
		treeCount = 0;
		Set<String> taxaToInclude = new LinkedHashSet<>();
		for (String t : prunedTree.getTaxaNames()) {
			taxaToInclude.add(t);
		}
		Set<String> before = new HashSet<>();
		collectTaxa(tree.getRoot(), before);
		Set<String> beforeMapped = new HashSet<>();
		for (String t : before) {
			beforeMapped.add(taxonMap.get(t));
		}
		before = beforeMapped;
		tree.setRoot(filter(tree.getRoot(), taxaToInclude, taxonMap));
		Set<String> after = new HashSet<>();
		collectTaxa(tree.getRoot(), after);
		before.removeAll(after);
		Log.warning("Removed: " + Arrays.toString(before.toArray()));
		
		TaxonSet ts = taxonSetForTree(tree.getRoot());
		CladeSet cladeSet1 = new CladeSet(tree, ts);
		
		
		while (srcTreeSet.hasNext()) {
			tree = srcTreeSet.next();
			if (tree == null) {
				return cladeSet1;
			}
			tree.setRoot(filter(tree.getRoot(), taxaToInclude, taxonMap));
			cladeSet1.add(tree);
			treeCount++;
		}
		return cladeSet1;
	}
	
	private void collectTaxa(Node node, Set<String> taxa) {
		if (node.isLeaf()) {
			taxa.add(node.getID());
		} else {
			for (Node child : node.getChildren()) {
				collectTaxa(child, taxa);
			}
		}
	}
	
	private TaxonSet taxonSetForTree(Node root) {
		List<Taxon> ts = new ArrayList<>();
		taxonSetForTree(tree.getRoot(), ts);

		TaxonSet tset = new TaxonSet();
		for (Taxon t : ts) {
			tset.taxonsetInput.get().add(t);
		}
		tset.initAndValidate();
		return tset;
	}

	private void taxonSetForTree(Node node, List<Taxon> ts) {
		if (node.isLeaf()) {
			ts.add(new Taxon(node.getID()));
		} else {
			for (Node child : node.getChildren()) {
				taxonSetForTree(child, ts);
			}
		}
	}

	/**
	 * rename leafs using taxonMap 
	 * filter out nodes from a binary tree if they are not in the taxaToInclude **/
	protected Node filter(Node node, Set<String> taxaToInclude, Map<String, String> taxonMap) {
		if (node.isLeaf()) {
			String id = node.getID();
			if (taxonMap.containsKey(id)) {
				node.setID(taxonMap.get(id));
			}
			if (taxaToInclude.contains(node.getID())) {
				
				return node;
			} else {
				return null;
			}
		} else {
			Node left_ = node.getLeft(); 
			Node right_ = node.getRight(); 
			left_ = filter(left_, taxaToInclude, taxonMap);
			right_ = filter(right_, taxaToInclude, taxonMap);
			if (left_ == null && right_ == null) {
				return null;
			}
			if (left_ == null) {
				return right_;
			}
			if (right_ == null) {
				return left_;
			}
			node.removeAllChildren(false);
			node.addChild(left_);
			node.addChild(right_);
			return node;
		}
	}	
	
	public static void main(String[] args) throws Exception {
		try {
			new Application(new GenerateLexicalConstraints(), "Generate Lexical Constraints", args);
		} catch (StackOverflowError e) {
			Log.warning("Are you running with enough stack space (e.g., try 'beast -Xss128m ...')?");
			throw e;
		} catch (Exception e) {
			throw e;
		}
	}

}
