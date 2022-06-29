#!/bin/sh
export BEAST_HOME=/Applications/BEAST\ 2.5.2/bin/
cd /Users/remco/data/beast/global/monos

"$BEAST_HOME/applauncher" GenerateLexicalConstraints -glottologTree dplace/iso.tree -xml dplace/geo-rc1197+almostnewwals+corrcals4.xml -treeID Tree.t:DPLACE -treeConfig treeset.cfg -burnin 10 -out /tmp/dplace.xml
