package monos;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import beast.app.beauti.BeautiDoc;
import beast.app.treeannotator.TreeAnnotator;
import beast.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beast.app.util.Application;
import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Input;
import beast.core.Logger;
import beast.core.Runnable;
import beast.core.Input.Validate;
import beast.core.util.CompoundDistribution;
import beast.core.util.Log;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.CladeSet;
import beast.evolution.tree.Tree;
import beast.math.distributions.AlmostMRCAPrior;
import beast.math.distributions.AlmostMultiMRCAPriors;
import beast.math.distributions.MRCAPrior;

@Description("test whether the trees in the treeConfig file used for GenerateLexicalConstraints can be parsed, "
		+ "without going through the complete GenerateLexicalConstraints process")
public class TestTreeFormat extends Runnable {

	/** inputs required for finding lexical analyses **/ 
	final public Input<File> treeConfigInput = new Input<>("treeConfig", 
			"configuration file that identifies lexical analyses. "
			+ "This should be a tab delimited file with three columns: name of the posterior tree set, "
			+ "name of the taxon map file, and scale factor for trees.", Validate.REQUIRED);

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		parseTrees();
	}

	private void parseTrees() throws IOException {
		
		String cfg = BeautiDoc.load(treeConfigInput.get());
		String [] strs = cfg.split("\n");
		for (String str : strs) {
			if (str.trim().length() > 0 && !str.trim().startsWith("#")) {
				String treefile = str.split("\t")[0];
				
				MemoryFriendlyTreeSet srcTreeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(treefile, 0);				
				srcTreeSet.reset();
				Tree tree = srcTreeSet.next();
				int treeCount = 0;
				System.err.print(treefile + " has... ");
				while (srcTreeSet.hasNext()) {
					try {
						tree = srcTreeSet.next();
					} catch (Throwable e) {
						e.printStackTrace();
						System.err.println("treeCount = " + treeCount);
						return;
					}
					treeCount++;
					if (treeCount % 100 == 0) { 
						System.err.print(".");
					}
//					if (treeCount == 999) {
//						System.err.print("|");
//					}
				}
				System.err.println(treeCount + " trees");
			}
		}		
	}
	
	
	public static void main(String[] args) throws Exception {
		new Application(new TestTreeFormat(), "Test Tree Format", args);

	}

}
