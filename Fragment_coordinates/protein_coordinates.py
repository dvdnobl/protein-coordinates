#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 22:04:52 2019

@author: davidnoble
"""

import pandas as pd
import pybedtools as bd

# accessing BED files of PacBio sequencing data and annotated genome
sequences = bd.BedTool('pacbio-190731-facs-assign.bed')
annotations = bd.BedTool('saccharomyces_cerevisiae.bed')

fragment_peaks = pd.read_csv('joint-frag-mle-peak.csv')

# finding overlap, where entirety of item from SEQUENCES has overlap with item in ANNOTATIONS
# some will still span introns, due to nature of ANNOTATIONS data
intersections = sequences.intersect(annotations, f=1.0, wo=True, nonamecheck=True).to_dataframe()

# creating DataFrames from BedTool objects (for ease)
sequences = sequences.to_dataframe(names = ['chrom', 'start', 'end', 
                                          'barcode', 'num_reads', 'strand'])
annotations = annotations.to_dataframe()

cols=[6,10,11,12,13,14]
intersects = intersections.drop(intersections.columns[cols], axis=1)
intersects.columns = ['chrom','start','end','barcode','num_reads',
                      'strand','geneStart','geneEnd','yorf','exons',
                      'exonLengths','exonStarts','overlap']

in_frame = intersects[(intersects['start'] - intersects['geneStart']) % 3 == 0]
print(in_frame.head()['exonStarts'])