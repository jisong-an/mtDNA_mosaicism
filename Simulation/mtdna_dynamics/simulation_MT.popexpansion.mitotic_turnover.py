###################################
#
#   Simulation : MT dynamics to infer mtDNA expansion (mitotic)
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
####################################



import numpy as np
import pandas as pd
import argparse
import sys
import random
from collections import Counter



def argument_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-g', '--generation', required=True, type=int, help='number of cell division')
	parser.add_argument('-m', '--mtcn', required=False, type=int, default=1000, help='mtDNA copy number')
	parser.add_argument('-n', '--ncells', required=False, type=int, default=1000, help='number of simulated cells')
	parser.add_argument('-i', '--iter', required=False, type=int, default=1, help='parameter sampling iteration')
	parser.add_argument('-o', '--output', required=True, help='output TSV file')
	args = vars(parser.parse_args())
	return args['generation'], args['mtcn'], args['ncells'], args['iter'], args['output']


def summarize_stat(cells, n_mitochondria):

	n_cells = cells.shape[0]
	freq_lst = []

	for k, lst in enumerate(cells, start=1):
		freq = np.array([count for _,count in Counter(lst).most_common(10)], dtype=np.float64)
		freq_uniq = np.array(sorted(list(set([count for _,count in Counter(lst).most_common(100)])), reverse=True)[:10], dtype=np.float64)
		freq_norm = freq/n_mitochondria
		freq_uniq_norm = freq_uniq/n_mitochondria

		freq_str = ':'.join(str(f) for f in freq_norm)
		freq_uniq_str = ':'.join(str(f) for f in freq_uniq_norm)
		freq_final = freq_str + '|' + freq_uniq_str
		freq_lst.append(freq_final)

	result = ['__'.join(str(f) for f in freq_lst)]

	return result


# define one turnover equivalent
def turnover_simul(cells,n_mitochondria, n_cells):

	for i in range(n_mitochondria):
		
		# Randomly select a column index to duplicate 
		duplicate_index = np.random.choice(n_mitochondria+i, size=(n_cells, 1), replace=True)
		#print(duplicate_index)
		rows = np.arange(n_cells)
		duplicated_column = cells[rows, duplicate_index[:,0]]
		cells = np.column_stack((cells, duplicated_column))

	# after all duplications, segregate randomly
	random_index = np.random.choice(2*n_mitochondria, size=n_mitochondria, replace=False)
	cells = cells[:, random_index]

	return cells


def cell_division(n_mitochondria, generation, n_cells, cells):

	# initial VAF
	result = summarize_stat(cells, n_mitochondria)
	data = {'generation' : [0], 'mtCN' : [f'{n_mitochondria}'], 'n_cells' : [f'{n_cells}'], 'result' : result}

	# cell division start
	for i in range(generation):

		# shuffle
		np.apply_along_axis(np.random.shuffle, axis=1, arr=cells)

		# cell turnover
		cells = turnover_simul(cells, n_mitochondria, n_cells)
		print(f'{i} th turnover end!')
		sys.stdout.flush()

		# record cell division info
		result = summarize_stat(cells, n_mitochondria)
		data['generation'].extend([f'{i+1}'])
		data['mtCN'].extend([f'{n_mitochondria}'])
		data['n_cells'].extend([f'{n_cells}'])
		data['result'].extend(result)
		#print(data)

	# result to dataframe
	df = pd.DataFrame(data)
	
	return df


def initVAF(n_mitochondria, n_cells):


	# make mutation to all mtDNA (1~n_mitochondria)
	mutant_temp = np.arange(1,n_mitochondria+1)
	cell_temp = np.tile(mutant_temp, n_cells)
	cell_vaf = np.reshape(cell_temp, (n_cells, n_mitochondria))

	# shuffle
	np.apply_along_axis(np.random.shuffle, axis=1, arr=cell_vaf)

	return cell_vaf



def main():

	generation, n_mitochondria, n_cells, iter, output = argument_parser()

	# simulate cell division
	for n in range(iter):

		# cells with many mutations
		cells_vaf = initVAF(n_mitochondria, n_cells)
		result = cell_division(n_mitochondria, generation, n_cells, cells_vaf)

		# split columns
		col_names = ['cell' + str(i) for i in range(1,n_cells+1)]
		result[col_names] = result['result'].str.split('__', expand=True)
		result = result.drop("result", axis=1)

		result.to_csv(output.replace('.tsv',f'.iter{n+1}.tsv'), sep='\t', index=False)


if __name__=='__main__':
	main()
	
