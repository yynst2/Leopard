from .misc import *

class motif(object):

	motif_path=''

	def __init__(self):
		self.length = 0
		self.matrix = []
		self.r_matrix = []
		self.score_matrix = []  # no longer useful as we'll put log transformed percentage score into the matrix
		self.tf = ''
		self.motif_id = ''
	def load_motif(self, motif_path):
		'''
		:param motif_path:
		:return:
		'''
		file = os.path.join(motif_path, self.motif_id + '.txt')
		df = pd.read_csv(file, header=0, index_col=0, sep='\t')
		self.length = df.shape[0]
		mat = df.values
		mat2 = copy.copy(mat)
		max_value = 0.0
		for i in range(df.shape[0]):
			max_value_temp = 0.0
			for j in range(df.shape[1]):
				if mat[i, j] < 0.001:
					mat[i, j] = 0.001
				mat2[i, j] = np.log2(mat[i, j] / 0.25)
				if mat2[i, j] >= max_value_temp:
					max_value_temp = mat2[i, j]
			max_value += max_value_temp
		self.matrix = mat2 / max_value
		print(self.motif_id, max_value, '\t')
		self.r_matrix = self.matrix.T[::-1].T[::-1]
	@staticmethod
	def load_all_motifs(motif_path, motif_mapping_file):
		'''
		:param motif_path:
		:param motif_dir:
		:param motif_mapping_file:    motif_id<T>TF<T>Index
		:return:
		'''
		motif_vec = []
		with open(file=motif_mapping_file, mode='r') as f:
			data = f.readlines()
		for _ in data:
			_ = _.rstrip('\n')
			m = motif()
			m.motif_id = _.split('\t')[0]
			m.tf = _.split('\t')[1]

			m.load_motif(motif_path=motif_path)
			motif_vec.append(m)
		return motif_vec
	@staticmethod
	def load_all_motifs_leopard(motif_path, motif_list, tf):
		'''
		:param motif_path:
		:param motif_dir:
		:param motif_mapping_file:    motif_id<T>TF<T>Index
		:return:
		'''
		motif_vec = []
		for _ in motif_list:
			m = motif()
			m.motif_id = _
			m.tf = tf

			m.load_motif(motif_path=motif_path)
			motif_vec.append(m)
		return motif_vec
	@staticmethod
	def standardize_motifs(motif_vec): # need test (todo)
		'''
		return zero padded motif vectors (including forward and reverse) for a list of known motifs
		'''

		# matrix size
		matrix_size=2*int(len(motif_vec))  # including forward and reverse
		# set max length
		max_length=0
		for i in motif_vec:
			if i.length>=max_length:
				max_length=i.length
		tmp_container=[]
		for i in motif_vec:
			# forward
			tmp_matrix=np.zeros(shape=(max_length,4))
			for x in range(i.length):
				for y in range(4):
					tmp_matrix[x,y]=i.matrix[x,y]
			tmp_matrix=tmp_matrix.reshape((max_length, 4, 1))
			tmp_container.append(tmp_matrix)
			# reverse
			tmp_matrix = np.zeros(shape=(max_length, 4))
			for x in range(i.length):
				for y in range(4):
					tmp_matrix[x, y] = i.r_matrix[x, y]
			tmp_matrix=tmp_matrix.reshape((max_length, 4, 1))
			tmp_container.append(tmp_matrix)
		# zero padding
		standardized_motifs=np.moveaxis(tmp_container,0,2)
		standardized_motifs= standardized_motifs.reshape((1,max_length,4,matrix_size))

		return standardized_motifs

	def __eq__(self, other):
		return other == self.motif_id
