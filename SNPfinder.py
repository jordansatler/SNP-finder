#!/usr/bin/env python

"""
Locates SNPs in locus files and concatenates into a
data set. Locus files assumed to be in nexus format.
Can retain all SNPs (linked) or can randomly grab a 
single SNP per locus (unlinked). Indels are not 
recognized as variable sites. Output is in phylip format.

date: 1 Oct 2018
author: J. Satler
version: 1

usage:
    python SNPfinder.py /path/to/nexus/files linked|unlinked
"""

import os
import sys
import random

def locus_files(files):
    """get list of nexus files"""
    return [os.path.join(files, f) for f in os.listdir(files)
            if f.endswith("nexus")]

def read_data(locus):
    """parse nexus file"""
    with open(locus, 'r') as l:
        start = -1
        m = {}
        for line in l:
            line = line.strip()
            if line == "matrix":
                start = 0
                continue
            elif line == ";":
                # return as list of lists for next function
                return [[k,''.join(v)] for (k,v) in m.items()]
            elif start == -1:
                continue
            else:
                if line:
                    if line.split()[0] in m:
                        # date are interleaved
                        m[line.split()[0]].extend(line.split()[1])
                    else:
                        # data are not interleaved
                        m[line.split()[0]] = [line.split()[1]]

def find_snps(mat):
    """locate snps in locus"""
    s = {n[0]:[] for n in mat}
    z_mat = [f[1] for f in mat]

    # find snps with allowed base pairs
    allowed = ['A', 'C', 'G', 'T']
    for index, i in enumerate(zip(*z_mat)):
        if len(set(i).intersection(allowed)) > 1:
            # found a snp
            for j in range(len(mat)):
                s[mat[j][0]].append(z_mat[j][index])
    return s

def get_random_snp(loc_snp):
    """select one snp at random"""
    unlinked_snp = random.randint(0, len(list(loc_snp.values())[0]) - 1)
    return {k:v[unlinked_snp] for (k,v) in loc_snp.items()}

def build_concat_matrix(taxa, snp_loci, linkage):
    "build concatenated snp matrix"
    snp_mat = {t:[] for t in taxa}
    for index, locus in sorted(snp_loci.items()):
        for i in taxa:
            if i in locus:
                snp_mat[i].extend(locus[i])
            else:
                missing = '?' * len(list(locus.values())[0])
                snp_mat[i].extend(missing)
    # write to file
    write_out(snp_mat, linkage)

def write_out(snp_mat, linkage):
    """write snp locus to file"""
    with open("snp_concat_" + linkage + ".phy", 'w') as out:
        out.write("{0} {1}\n".format(len(snp_mat),
                                    len(list(snp_mat.values())[0])))
        for k, v in sorted(snp_mat.items()):
            out.write("{0}{1}{2}\n".format(k, ' ' * (30 - len(k))
                                            ,''.join(v)))

def main():
    if len(sys.argv) != 3:
        print("python SNPfinder.py /path/to/nexus/files linked|unlinked")
        sys.exit()

    snp_loci = {}
    linkage = "linked"

    data = locus_files(sys.argv[1])
    for index, i in enumerate(data):
        loc = read_data(i)
        loc_snps = find_snps(loc)

        # check if locus has any snps
        if not list(loc_snps.values())[0]:
            continue

        # if using unlinked snps
        if sys.argv[2].lower().startswith("u"):
            linkage = "unlinked"
            loc_snps_unlinked = get_random_snp(loc_snps)
            # add to master dictionary
            snp_loci[index] = loc_snps_unlinked
            continue

        # add linked snp matrix to dictionary if using linked snps
        snp_loci[index] = loc_snps

    # get full taxon names
    taxa = []
    for k, v in snp_loci.items():
        for name in v.keys():
            if name not in taxa:
                taxa.append(name)

    # build and write concatenated matrix
    snp_mat = build_concat_matrix(taxa, snp_loci, linkage)

if __name__ == '__main__':
    main()