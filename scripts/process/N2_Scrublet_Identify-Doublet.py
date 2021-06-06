#!/usr/bin/env python

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import getopt


def main(argv):
   inputdir = ''
   edr = ''
   try:
      opts, args = getopt.getopt(argv,"hi:r:",["inputdir=","edr="])
   except getopt.GetoptError:
      print('test.py -i <inputdir> -r <expected doublet ratio>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('test.py -i <inputdir> -r <expected doublet ratio>')
         sys.exit()
      elif opt in ("-i", "--inputdir"):
         inputdir = arg
      elif opt in ("-r", "--edr"):
         edr = arg

# change the dir to your input data directory 
   os.chdir('/lustre1/lch3000_pkuhpc/liuyt/SP/D200914-scRNA-data/cellranger-count/'+inputdir)
   sample = os.path.basename(os.getcwd())

 #  plt.rcParams['font.family'] = 'sans-serif'
   plt.rcParams['font.sans-serif'] = 'Arial'
   plt.rc('font', size=14)
   plt.rcParams['pdf.fonttype'] = 42

   input_dir = 'filtered_feature_bc_matrix'
   counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
   genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1))
   print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
   print('Number of genes in gene list: {}'.format(len(genes)))

   print('EDR: {}'.format(edr))
   scrub = scr.Scrublet(counts_matrix, expected_doublet_rate = float(edr))

   doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                          min_cells=3,
                                                          min_gene_variability_pctl=85,
                                                          n_prin_comps=30)

   scrub.plot_histogram();
   plt.savefig(sample + '_scrublet_EDR' + edr + '_hist.pdf')

   print('Running UMAP...')
   scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
   scrub.plot_embedding('UMAP', order_points=True)
   plt.savefig(sample + '_scrublet_EDR0.02_umap.pdf')

   outf1 = sample + "_scrublet_EDR" + edr + "_DoubletScores.csv"
   outf2 = sample + "_scrublet_EDR" + edr + "_PredictedDoublets.csv"
   np.savetxt(outf1, doublet_scores, delimiter=',')
   np.savetxt(outf2, predicted_doublets, delimiter=',')


if __name__ == "__main__":
   main(sys.argv[1:])
