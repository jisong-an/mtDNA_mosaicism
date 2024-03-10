###################################
#
#   Simulation : MT dynamics to infer mtDNA heteroplasmy expansion (mitotic)
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
	split_cells = [list(filter(None,sublist)) for sublist in np.char.split(cells,':')]
	freq_lst = []

	for k, lst in enumerate(split_cells, start=1):
		var = np.array([int(i) for i in np.concatenate(lst) if i != ''], dtype=np.int64)
		freq = np.array(list(Counter(var).values()), dtype=np.float64)
		freq_norm = freq/n_mitochondria
		freq_norm_sort = sorted(freq_norm,reverse=True)[0:9]

		freq_str = ':'.join(str(f) for f in freq_norm_sort)
		freq_uniq_str = ':'.join(str(f) for f in freq_norm_sort)
		freq_final = freq_str + '|' + freq_uniq_str
		freq_lst.append(freq_final)


	result = ['__'.join(str(f) for f in freq_lst)]

	return result



def cell_division(n_mitochondria, generation, n_cells, mrate_gen):

	mut_idx = 0
	cells = np.zeros((n_cells, n_mitochondria*2), dtype='<U21')
	turnover = generation

	# cell division start
	for i in range(turnover):

		for j in range(n_mitochondria):

			# Generate mutations
			mut = np.random.poisson(mrate_gen, (n_cells,))
			mut_random = np.random.choice(n_cells,1)
			mut_common = mut[mut_random][0]
			
			duplicate_index = np.random.randint(n_mitochondria+j)

			# replace mtDNA & insert mutation (mutation occurs during duplication)
			if mut_common > 0:
				rows = np.arange(n_cells)
				cols = np.array([n_mitochondria+j])  # newly synthesized mtDNA
				cells[:,n_mitochondria+j] = cells[:,duplicate_index]
				values = np.core.defchararray.add(cells[rows,cols], np.array([f':{mut_idx}'], dtype='<U21'))  # add string to specific elements
				np.put(cells, np.ravel_multi_index((rows[:, None], cols[:, None]), cells.shape), values)  # update cells (put(array,indices, values, mode='raise'))
				mut_idx += 1

			# when no mutation occurs, only turnover
			if mut_common == 0:
				cells[:,n_mitochondria+j] = cells[:,duplicate_index]
		

		# shuffle when one turnover ends
		np.apply_along_axis(np.random.shuffle, axis=1, arr=cells)

		random_index = np.random.choice(2*n_mitochondria, size=2*n_mitochondria, replace=False)
		cells = cells[:, random_index]	


		# record cell division info
		if (i+1)%10==0:

			select_cells = cells[:,0:n_mitochondria]
			
			# select clones & summary statistics
			result = summarize_stat(select_cells, n_mitochondria)
			print(i+1)
			sys.stdout.flush()
			
			if i+1==10:
				data = {'generation' : [f'{i+1}'], 'mtCN' : [f'{n_mitochondria}'], 'n_cells' : [f'{n_cells}'], 'result' : result}

			else:
				data['generation'].extend([f'{i+1}'])
				data['mtCN'].extend([f'{n_mitochondria}'])
				data['n_cells'].extend([f'{n_cells}'])
				data['result'].extend(result)

	# result to dataframe
	df = pd.DataFrame(data)
	

	return df




def main():

	generation, n_mitochondria, n_cells, iter, output = argument_parser()

	mrate     = 5e-08
	mrate_gen = 16569*mrate
	print(mrate_gen)

	# simulate cell division
	for n in range(iter):

		# cells with many mutations
		result = cell_division(n_mitochondria, generation, n_cells, mrate_gen)

		# split columns
		col_names = ['cell' + str(i) for i in range(1,n_cells+1)]
		result[col_names] = result['result'].str.split('__', expand=True)
		result = result.drop("result", axis=1)

		result.to_csv(output.replace('.tsv',f'.iter{n+1}.tsv'), sep='\t', index=False)


if __name__=='__main__':
	main()
	
