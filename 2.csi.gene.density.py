#!/usr/bin/python
# coding=utf-8
# Created by yinzhaoping on 2017.09.07
# get chr density


def read_gene(genefile):
    s = []
    chrms = []
    with open(genefile) as lines:
        for line in lines:
            if line.find('gene') >= 0:
                t = line.strip().split('\t')
                if t[0].find('Un') >= 0:
                    continue
                chrom = int(t[0].split('hr')[1])
                start = int(t[3])
                end = int(t[4])
                gene = t[8].split(';')[0].replace('ID=', '')
                b = [chrom, gene, start, end]
                s.append(b)
                chrms.append(chrom)
    s.sort(key=lambda x: (x[0], x[2]))
    m = list(set(chrms))
    m.sort()
    return s, m


genefile = 'csi.gene.models.gff3'
outfile = 'csi_gene_density.txt'
f = open(outfile, 'w')
section = 100
domain = 200000
genes, chrms = read_gene(genefile)

start = 0
for i, ch in enumerate(chrms):
    if ch > 9:
        continue
    print(ch)

    ch_tag = 0
    for j in range(0, 1000):
        a = domain*j
        b = a + domain - 1
        s_tag = 0
        n = 0
        while 1:
            for k in range(start, start+section):
                if k >= len(genes):
                    s_tag = 1
                    break
                if genes[k][2] > b:
                    s_tag = 1
                    break
                if genes[k][0] != ch:
                    # print start, genes[k][0], ch
                    ch_tag = 1
                    break
                if a <= genes[k][2] <= b:
                    start = k+1
                    length = genes[k][3] - genes[k][2] + 1
                    # n += 1
                    n += length
                    # print ch, a, b, genes[k],start
            if s_tag or ch_tag:
                break

        f.write('chr%d %d %d %d\n' % (ch, a, b, n))
        if ch_tag == 1:
            break
f.close()
