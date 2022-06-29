This folder contains data for creating lexical constraints using the cd GenerateLexicalConstraints app

You can run the GenerateLexicalConstraints app from this directory using (see README.md in /monos for more details)
applauncher GenerateLexicalConstraints -glottologTree all-glotto.tree -xml ../../output/edge.xml -treeID Tree.t:EDGE -treeConfig treeset-glotto-hand-corrected.cfg -burnin 10 -out /tmp/edge-with-lexical.xml
This assumes that ../../output/edge.xml is the BEAST XML file for which lexical constraints need to be generated.


Files:
./treeset-glotto-hand-corrected.cfg configuration file used by GenerateLexicalConstraints. This should be a tab delimited file with three columns: name of the posterior tree set, name of the taxon map file, and scale factor for trees


./d-atlantic-congo/glotto-corrected.ac.txt
./d-austronesian/glotto-corrected.au.txt
./d-indo-european/glotto-corrected.ie.txt
./d-turkic/glotto.turkic.txt
./d-uralic/glotto-corrected.uralic.txt
./dravidian/glotto.drav.txt
./paman-nyungan/glotto-corrected.pn.txt
./semitic/glotto-corrected.semitic.txt
./sino-tibetan/glotto-corrected.st.txt

The above files map glotto codes to the taxon names as used in the tree files.


./d-atlantic-congo/posterior.trees
./d-austronesian/posterior.trees
./d-indo-european/posterior.trees
./d-turkic/posterior.trees
./d-uralic/posterior.trees
./dravidian/posterior.trees
./paman-nyungan/posterior.trees
./semitic/posterior.trees
./sino-tibetan/posterior.trees

The above files are sourced from D-PLACE github repop https://github.com/D-PLACE/dplace-data/


all-glotto.tree glottolog tree (in monos/data/tree_glottolog_newick.txt) where all meta data is removed using search/replace in a text editor. The original tree puts families trees on one per line, this tree replaces all newlines with commas, and wraps the result in brackets to create a valid tree.



