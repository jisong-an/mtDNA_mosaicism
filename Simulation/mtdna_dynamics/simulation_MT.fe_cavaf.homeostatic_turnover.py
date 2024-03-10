###################################
#
#   Simulation : MT dynamics with initial VAF (homeostatic)
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
	parser.add_argument('-g', '--generation', required=True, type=int, help='number of turnover')
	parser.add_argument('-v', '--initvaf', required=True, type=float, help='caVAF of fertilized-egg origin variant')
	parser.add_argument('-m', '--mtcn', required=False, type=int, default=1000, help='mtDNA copy number')
	parser.add_argument('-n', '--ncells', required=False, type=int, default=1000, help='number of simulated cells')
	parser.add_argument('-i', '--iter', required=False, type=int, default=1, help='parameter sampling iteration')
	parser.add_argument('-o', '--output', required=True, help='output TSV file')
	args = vars(parser.parse_args())
	return args['generation'], args['initvaf'], args['mtcn'], args['ncells'], args['iter'], args['output']


def summarize_stat(cells, n_mitochondria):

	n_cells = cells.shape[0]
	freq_lst = []

	for k, lst in enumerate(cells, start=1):

		freq = round(sum(lst)/n_mitochondria,3)
		freq_lst.append(freq)

	# wild-type & heteroplasmy & homoplasmy
	wt = round(sum(x < 0.005 for x in freq_lst)/n_cells,3)
	hetero = round(sum(0.005 <= x < 0.9 for x in freq_lst)/n_cells,3)
	homo = round(sum(x >= 0.9 for x in freq_lst)/n_cells,3)

	# mean & sd
	mean_val = round(np.mean(freq_lst),3)
	std_val = round(np.std(freq_lst),3)

	# VAF count in specific range
	bin_counts = np.histogram(freq_lst,bins=np.concatenate((np.array([0.005, 0.05]),np.arange(0.1, 1.05, 0.05))))[0]
	bin_counts_str = ':'.join(str(count) for count in bin_counts)

	freq_lst.extend([f'{wt}',f'{hetero}',f'{homo}', f'{mean_val}', f'{std_val}', f'{bin_counts_str}'])
	result = ['__'.join(str(f) for f in freq_lst)]

	return result


# define one turnover equivalent
def turnover_simul(cells,n_mitochondria):

	num_columns = cells.shape[1]
	
	for _ in range(n_mitochondria):
		
		# Randomly select a column index to duplicate and remove
		duplicate_index = np.random.randint(num_columns)
		remove_index = np.random.randint(num_columns)
		cells[:, remove_index] = cells[:, duplicate_index]
		
	return cells



def cell_division(n_mitochondria, generation, initvaf, n_cells, cells):


	# initial VAF
	result = summarize_stat(cells, n_mitochondria)
	data = {'generation' : [0], 'initVAF' : [f'{initvaf}'], 'mtCN' : [f'{n_mitochondria}'], 'n_cells' : [f'{n_cells}'], 'result' : result}

	# for fixation
	fixed_data = {'cell' : [], 'fixVAF' : [], 'fixGen' : []}


	# cell division start
	for i in range(generation):

		# shuffle
		np.apply_along_axis(np.random.shuffle, axis=1, arr=cells)

		# cell turnover
		cells = turnover_simul(cells, n_mitochondria)
		print(f'{i} th turnover end!')
		sys.stdout.flush()

		# check fixation
		for j,cell in enumerate(cells):
			if np.all(cell==0) and j not in fixed_data['cell']:
				fixed_data['cell'].append(j)
				fixed_data['fixVAF'].append(0)
				fixed_data['fixGen'].append(i+1)
			elif np.all(cell==1) and j not in fixed_data['cell']:				
				fixed_data['cell'].append(j)
				fixed_data['fixVAF'].append(1)
				fixed_data['fixGen'].append(i+1)


		# record cell division info
		# select clones & summary statistics
		result = summarize_stat(cells, n_mitochondria)
		data['generation'].extend([f'{i+1}'])
		data['initVAF'].extend([f'{initvaf}'])
		data['mtCN'].extend([f'{n_mitochondria}'])
		data['n_cells'].extend([f'{n_cells}'])
		data['result'].extend(result)

	# result to dataframe
	df = pd.DataFrame(data)
	df2 = pd.DataFrame(fixed_data)
	

	return df,df2


def initVAF(n_mitochondria, initvaf, n_cells):


	# make mutant mtDNA
	n_mutant = round(initvaf*n_mitochondria)
	mutant_temp = np.concatenate((np.repeat(0,n_mitochondria-n_mutant), np.repeat(1,n_mutant)), axis=None)  # 0 : wildtype, 1 : mutant
	cell_temp = np.tile(mutant_temp, n_cells)
	cell_vaf = np.reshape(cell_temp, (n_cells, n_mitochondria))

	# shuffle
	np.apply_along_axis(np.random.shuffle, axis=1, arr=cell_vaf)

	return cell_vaf



def main():

	generation, initvaf, n_mitochondria, n_cells, iter, output = argument_parser()

	# simulate cell division
	for n in range(iter):

		print(f'{n} th iteration start!')
		sys.stdout.flush()

		# cells after bottleneck
		cells_vaf = initVAF(n_mitochondria, initvaf, n_cells)
		result, fixed_result = cell_division(n_mitochondria, generation, initvaf, n_cells, cells_vaf)

		# split columns
		col_names = ['cell' + str(i) for i in range(1,n_cells+1)]
		col_names += ['wt_ratio','hetero_ratio','homo_ratio','meanVAF', 'sdVAF', 'VAFdist']
		result[col_names] = result['result'].str.split('__', expand=True)
		result = result.drop("result", axis=1)

		result.to_csv(output.replace('.tsv',f'.iter{n+1}.tsv'), sep='\t', index=False)
		fixed_result.to_csv(output.replace('.tsv',f'.fixdata.iter{n+1}.tsv'), sep='\t', index=False)


if __name__=='__main__':
	main()
	
