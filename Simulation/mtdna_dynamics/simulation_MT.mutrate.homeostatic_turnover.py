###################################
#
#   Simulation : MT dynamics -estimate mutation rate (homeostatic)
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
####################################



import numpy as np
import pandas as pd
import argparse
import sys
from collections import Counter


def argument_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-t', '--turnover', required=True, type=float, help='turnover rate per year')
	parser.add_argument('-a', '--age', required=True, type=float, help='age')
	parser.add_argument('-c', '--mtcn', required=True, type=int, help='mtDNA copy number')
	parser.add_argument('-s', '--sampling', required=True, type=int, help='number of clones to sample')
	parser.add_argument('-m', '--maxcell', required=False, type=int, default=1000, help='maximum cells (limit)')
	parser.add_argument('-n', '--sampleiter', required=False, type=int, default=1, help='sampling iteration')
	parser.add_argument('-i', '--totaliter', required=False, type=int, default=1, help='total iteration')
	parser.add_argument('-o', '--output', required=True, help='output TSV file')
	args = vars(parser.parse_args())
	return args['turnover'], args['age'], args['mtcn'], args['sampling'], args['maxcell'], args['sampleiter'], args['totaliter'], args['output']


def simulation_sequencing(cells, sampling, n_mitochondria):

	n_cells = cells.shape[0]
	sample_indice = np.random.choice(range(n_cells), size=sampling, replace=False)
	cells_sampling = cells[sample_indice,:]


	# VAF distribution (summary statistics)
	split_cells = [list(filter(None,sublist)) for sublist in np.char.split(cells_sampling,':')]
	freq_lst = []

	for k, lst in enumerate(split_cells, start=1):
		var = np.array([int(i) for i in np.concatenate(lst) if i != ''], dtype=np.int64)
		freq = np.array(list(Counter(var).values()), dtype=np.float64)
		freq_norm = freq/n_mitochondria
		
		freq_temp = []
		for i in np.arange(0,1,0.05):
			if i==0: 
				i=0.005
				j=0.05
			else:
				j=i+0.05
			count = np.count_nonzero((freq_norm > i) & (freq_norm <= j)) # variant count in specific arrange
			freq_temp.append(str(count))
		freq_temp.append(str(np.count_nonzero(freq_norm > 0.005))) # total variant count (>0.5%)
		freq_lst.append(f'sample{k}|' + ':'.join(freq_temp))

	result = ['__'.join(freq_lst)]


	return result



def cell_division(n_mitochondria, generation, max_cell, mrate, mrate_gen, sampling, sampleiter, output):

	n_cells = max_cell
	mut_idx = 0
	cells = np.zeros((n_cells, n_mitochondria), dtype='<U21')
	turnover = generation*n_mitochondria  # total replication->degradation times

	for i in range(turnover):

		if (i+1)%n_mitochondria==0:  # one turnover equivalent ends, shuffle MT (if shuffle per turnover, it takes long time)
			# shuffle
			np.apply_along_axis(np.random.shuffle, axis=1, arr=cells)

		# Generate mutations
		mut = np.random.poisson(mrate_gen, (n_cells,))
		mut_random = np.random.choice(n_cells,1)
		mut_common = mut[mut_random][0]

		# random choose turnover index 
		duplicate_index = np.random.choice(n_mitochondria, size=(n_cells, 1), replace=True)
		remove_index = np.random.choice(n_mitochondria, size=(n_cells, 1), replace=True)

		# replace mtDNA & insert mutation (mutation occurs during duplication)
		if mut_common > 0:
			rows = np.arange(n_cells)
			cells[rows, remove_index[:,0]] = cells[rows, duplicate_index[:,0]]  # one mtDNA duplicates & one mtDNA degrades
			cols = remove_index[:,0]  # newly synthesized mtDNA
			values = np.core.defchararray.add(cells[rows,cols], np.array([f':{mut_idx}'], dtype='<U21'))  # add string to specific elements
			np.put(cells, np.ravel_multi_index((rows[:, None], cols[:, None]), cells.shape), values)  # update cells (put(array,indices, values, mode='raise'))
			mut_idx += 1

		# when no mutation occurs, only turnover
		if mut_common == 0:
			rows = np.arange(n_cells)
			cells[rows, remove_index[:,0]] = cells[rows, duplicate_index[:,0]]  # one mtDNA duplicates & one mtDNA degrades

	return cells



def main():

	turnover, age, mtcn, sampling, max_cell, sampleiter, totaliter, output = argument_parser()
	n_mitochondria = mtcn
	generation = round(turnover*age)

	for n in range(totaliter):

		print(f'{n} th iteration start!')
		
		# select log value first and then change
		mrate_log = np.random.uniform(-9,-3)  # 10^-9 to 10^-3
		mrate     = pow(10,mrate_log)
		mrate_gen = 16569*mrate
		
		print(f'    mrate : {mrate}')
		print(f'    mrate_gen : {mrate_gen}')
		sys.stdout.flush()

		# simulate cell division
		cells = cell_division(n_mitochondria, generation, max_cell, mrate, mrate_gen, sampling, sampleiter, output)	

		# simulate sequencing
		for m in range(sampleiter):
			
			result = simulation_sequencing(cells, sampling, n_mitochondria)

			if m==0:
				data = {'turnover' : [f'{turnover}'], 'age' : [f'{age}'], 'mtCN' : [f'{mtcn}'], 'generation' : [f'{generation}'], 'sampling': [f'{sampling}'], 'max_cell' : [f'{max_cell}'], 'mrate' : [f'{mrate}'], 'mrate_gen' : [f'{mrate_gen}'], 'result' : result}
			else:
				data['turnover'].extend([f'{turnover}'])
				data['age'].extend([f'{age}'])
				data['mtCN'].extend([f'{mtcn}'])
				data['generation'].extend([f'{generation}'])
				data['sampling'].extend([f'{sampling}'])
				data['max_cell'].extend([f'{max_cell}'])
				data['mrate'].extend([f'{mrate}'])
				data['mrate_gen'].extend([f'{mrate_gen}'])
				data['result'].extend(result)

		# result to dataframe
		df = pd.DataFrame(data)
		df.to_csv(output.replace('.tsv',f'.iter{n+1}.tsv'), sep='\t', index=False)


if __name__=='__main__':
	main()
	
