# Skipping project
Genome assembly problem is like skipping. Many many rounds, again and again.
# Current status
* Build kmer count as vertices (31-mer) (on ecoli) (&#8730;)
* Filter singleton kmer (&#8730;)
* Glue non-branching kmer path (&#8730;)
* Build (small k)-mer graph to early filter (big k)-mer (on fruit fly and even human) (&#8730;)
* Simplify the non-branching path graph (&#8730;)
* What is the correct "weight" of each vertice on graph? (&#8730;)
* Extend non-branching path through simple erroneous branching (&#8730;)
* Resolve simple bubble (bubble, yes bubble) <----- current
* Check some miss-assembly contigs even when only simply non-branching path

...............................

# Side toys project
* Get subgraph around a node within some distances
