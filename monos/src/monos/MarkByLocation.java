package monos;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.PrintStream;
import java.util.*;

import javax.imageio.ImageIO;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Runnable;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import babel.tools.Nexus2Newick;

@Description("Mark tree by internal node location")
public class MarkByLocation extends Runnable {
	final public Input<TreeFile> treesInput = new Input<>("trees", "NEXUS file containing a tree set",
			Validate.REQUIRED);
	final public Input<OutFile> outputInput = new Input<>("out", "output file, or stdout if not specified",
			new OutFile("[[none]]"));
	final public Input<List<File>> locationFileInput = new Input<>("locationFile", "specify locations in format "
			+ "$(latitude) $(longitude) one per line",
			new ArrayList<>());
	final public Input<String> tagInput = new Input<>("tag", "meta-data tag used to mark clade",
			"clade");


	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		int width = 1024;
		int height = 1024;
		width = 512;
		height = 512;
//		width = 100;
//		height = 100;
		BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		Graphics g = image.getGraphics();
	
		Log.warning("Generating dot map");
		int n = locationFileInput.get().size();
		int i = 0;
		int [] colours = new int[n];
		// map colour to file name
		Map<Integer,String> map = new HashMap<>();
		for (File f : locationFileInput.get()) {
			int c = Color.HSBtoRGB((i + 0.0f)/n, 0.9f, 0.9f);
			map.put(c, f.getName().substring(0,f.getName().indexOf('.')));
			colours[i] = c;
			g.setColor(new Color(c));
			String str = BeautiDoc.load(f);
			String [] strs = str.split("\n");
			for (String s : strs) {
				String [] strs2 = s.split("\\s+");
				double latitude = Double.parseDouble(strs2[0]);
				double longitude = Double.parseDouble(strs2[1]);
				int x = (int)((180+longitude) * width / 360.0);
				int y = height - (int)((90+latitude) * height / 180.0);
				g.fillOval(x-1, y-1, 3, 3);
			}
			i++;
		}
		ImageIO.write(image, "png", new File("/tmp/MarkByLocation.png"));

		Log.warning("Generating closest region map");
		int [] rgb = new int[width * height];
		image.getRGB(0, 0, width, height, rgb, 0, width);
		for (int x = 0; x < width; x++) {
			for (int y = 0; y < height; y++) {
				int c = rgb[x + y*width];
				if (c == 0) {
					 int x0 = 0, y0 = 0;
					 int dx = 0;
					 int dy = -1;
					 do {
				        if (x0 == y0 || (x0 < 0 && x0 == -y0) || (x0 > 0 && x0 == 1-y0)) {
						  	int tmp = dx;
						   	dx = -dy;
						   	dy = tmp;
						}
						x0 = x0+dx;
						y0 = y0+dy;
					 } while (x+x0 < 0 || x+x0 >= width || y+y0 < 0 || y+y0 >= height || rgb[x+x0+(y+y0)*width] == 0);
					 g.setColor(new Color(rgb[x+x0+(y+y0)*width]));
					 g.fillOval(x-1, y-1, 3, 3);
				}
			}
			
			if (x*1000/width %10== 0) {System.err.print(".");}
		}
		System.err.println(".");
		ImageIO.write(image, "png", new File("/tmp/MarkByLocation2.png"));
		image.getRGB(0, 0, width, height, rgb, 0, width);
		
		// open file for writing
		PrintStream out = System.out;
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			out = new PrintStream(outputInput.get());
			Log.warning("Writing to file " + outputInput.get().getPath());
		}

		// read trees one by one, relabel and write out relabeled tree in newick format
		MemoryFriendlyTreeSet trees = new TreeAnnotator().new MemoryFriendlyTreeSet(treesInput.get().getAbsolutePath(),
				0);
		trees.reset();
		Tree tree = trees.next();
		trees.reset();
		String tag = tagInput.get();
		int k = 0;
		while (trees.hasNext()) {
			tree = trees.next();
			for (Node node : tree.getNodesAsArray()) {
				Object o = node.getMetaData("location");
				if (o != null && o instanceof Double[]) {
					Double [] location = (Double[]) o;
					double longitude = location[1];
					double latitude = location[0];
					int x = (int)((180+longitude) * width / 360.0);
					int y = height - (int)((90+latitude) * height / 180.0);
					if (y == height) {y--;};
					if (x+y*width >= rgb.length) {
						k++;
						x=rgb.length-1;y=0;
					}
					int colour = rgb[x+y*width];
					node.metaDataString=tag + "=\"" + map.get(colour) + "\""; 
				}
			}
			StringBuilder buf = new StringBuilder();
			Nexus2Newick.toShortNewick(tree.getRoot(), buf, true);
			out.println(buf.toString());
		}
		out.println();
		Log.warning(k + " corrections");

		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			out.close();
		}
		
		Log.warning("Done");
	}

	public static void main(String[] args) throws Exception {
		new Application(new MarkByLocation(), "Mark by location", args);
	}

}
