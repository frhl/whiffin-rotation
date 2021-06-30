#!/usr/bin/env python3


# from https://link.springer.com/article/10.1186/1471-2105-9-192
from ushuffle import shuffle, Shuffler


def shuffle_utrs(seq, k = 2, max_iter = 10):
	
	seq = seq.encode('UTF-8')
	#seq = b"ababcdcdabdcabvababab"
	shuffler = Shuffler(seq, k)
	mat = list()
	for i in range(max_iter):
	    seqres = shuffler.shuffle()
	    seqres = seqres.decode('UTF-8')
	    mat.append(seqres)
	return(mat)

	#seq = seq.encode('UTF-8')
	#shuffler = Shuffler(seq, 2)
	#shuffler.shuffle()
	#return(Shuffler(seq, 2))


