# SNP-finder

This script will concatenate SNPs from a set of loci. Loci are assumed to be in nexus format with ".nexus" file ending. All SNPs can be retained (linked) or a single SNP randomly selected per locus (unlinked) can be retained. Indels are not recognized as variable sites. Output file will be in phylip format.


usage:  
```python
    python SNPfinder.py /path/to/nexus/files linked|unlinked
```
