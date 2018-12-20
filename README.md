# Skipping project
Genome assembly problem is like skipping. Many many rounds, again and again.
# Current status
* Build kmer count as vertices (31-mer) (on ecoli) (&#8730;)
* Filter singleton kmer (&#8730;)
* Glue non-branching kmer path (&#8730;)
* Build (small k)-mer graph to early filter (big k)-mer (on fruit fly and even human) (&#8730;)
* Simplify the non-branching path graph <--------------------- current
---------> Resolve slow couting edges
---------> Resolve memory redundant
---------> Fix weird non-branching path not glued together
* Extend non-branching path through simple erroneous branching
* Resolve simple bubble (bubble, yes bubble)
...............................

# Side toys project
* Get subgraph around a node within some distances
