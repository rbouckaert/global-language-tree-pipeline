package monos;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beastfx.app.tools.Application;
import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Logger;
import beast.base.inference.Runnable;
import beast.base.core.Input.Validate;
import beast.base.inference.CompoundDistribution;
import beast.base.core.Log;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.CladeSet;
import beast.base.evolution.tree.Tree;
import almostbeast.math.distributions.AlmostMRCAPrior;
import almostbeast.math.distributions.AlmostMultiMRCAPriors;
import beast.base.evolution.tree.MRCAPrior;

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
