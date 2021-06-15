from numpy import zeros, array

def initialise_score_tb(n_x, n_y):
    """returns zeroed score and traceback matrices"""
    # recalling that the first index is "rows", the y-axis
    F = zeros((n_y + 1, n_x + 1), dtype=int)
    for x in range(1, n_x + 1):
        F[0, x] = -x

    for y in range(1, n_y + 1):
        F[y, 0] = -y

    # TB is an "object" matrix. I'm using it to store
    # coordinates for the best path entering each element
    TB = zeros(F.shape, dtype='O')
    # different initialision
    for y in range(TB.shape[0]):
        for x in range(TB.shape[1]):
            if y == 0 and x != 0:
                TB[y, x] = (0, x - 1)
            elif y != 0 and x == 0:
                TB[y, x] = (y - 1, 0)
            else:
                TB[y, x] = (0, 0)

    return F, TB

def nw_score_matrix(F, TB, X, Y, n_x, n_y):
    """returns the NW scored and traceback matrices"""
    for x in range(1, n_x + 1):
        base_1 = X[x - 1] # adjust for +1 size of M
        for y in range(1, n_y + 1):
            base_2 = Y[y - 1]
            score = 1 if base_1 == base_2 else -1
            match = F[y - 1, x - 1] + score, (y - 1, x - 1)
            delete = F[y - 1, x] - 1, (y - 1, x)
            insert = F[y, x - 1] - 1, (y, x - 1)
            F[y, x], TB[y, x] = max([match, insert, delete])
    return F, TB

def global_alignment(TB, X, Y):
    """returns Needleman-Wunsch aligned X and Y"""
    s1, s2 = [], []
    # the ending coordinate, recalling TB
    # is +1 on each dimension
    y, x = len(Y), len(X)
    while (y, x) != (0, 0):
        # we get the next coord in the traceback
        y_1, x_1 = TB[y, x]
        b1 = X[x - 1] if x != x_1 else '-'
        b2 = Y[y - 1] if y != y_1 else '-'
        y, x = y_1, x_1

        s1.append(b1)
        s2.append(b2)

    # as the alignment is built backwards, we reverse
    s1.reverse()
    s2.reverse()
    return "".join(s1), "".join(s2)


# Run the subroutines
def align(seq1,seq2,query_label="",subject_label="", verbose = False):
    ''' A function that runs NW-algorithm, and collate the 
    gathers some statistics and subequently prints it
    
    returns tuble with (length, identities, gaps)
    
    nicely '''
    
    X = seq1; Y = seq2

    # Run algorithm
    n_x, n_y = len(X), len(Y)
    F, TB = initialise_score_tb(n_x, n_y)
    F, TB = nw_score_matrix(F, TB, X, Y, n_x, n_y)
    
    # Produce some alignment stats
    aligned = global_alignment(TB, X, Y)
    X = aligned[0]; Y = aligned[1]
    
    # Produce some statistics
    Z = str()
    
    # Generate stats and specially, alignment metadata (I know this is big O intensive, but it looks neat.)
    align_match = 0; align_gap = 0
    for i in range(len(X)):
        if X[i] == Y[i]: 
            Z += "|"
            align_match += 1
        elif X[i] == "-" or Y[i] == "-":
            Z += " "
            align_gap += 1
        elif X[i] != Y[i]:
            Z += "*"
    
    total_nt = len(Z)
    
    if verbose:
        # Show some information about the alignment
        print("# Query = {0} ~ {1} nt".format(query_label,len(X.replace("-",""))))
        print("# Subject = {0} ~ {1} nt".format(subject_label,len(Y.replace("-",""))))
        print("# Alignment = {0} nt\n".format(len(aligned[0])))

        # Show matches/gaps and method
        info_match = " identities={0}/{1} ({2}%)".format(align_match,total_nt,(round(100*align_match/total_nt,2)))
        info_gap = " gaps={0}/{1} ({2}%)\n".format(align_gap,total_nt,(round(100*align_gap/total_nt,2)))
        print(" method: Needleman-wunch")
        print("{0}, {1}".format(info_match,info_gap))

        # Iterate through the alignment for a nice print
        for i in range(0,len(aligned[0]),60):

            # Generate the apropiate sequence
            xslice = len(X[i:i+60].replace("-","")); yslice = len(Y[i:i+60].replace("-",""))
            print("Query:\t {0} \t {1}".format(X[i:i+60],xslice))
            print("\t {0} \t".format(Z[i:i+60]))
            print("Sbjct:\t {0} \t {1}\n".format(Y[i:i+60],yslice))

    return total_nt, align_match, align_gap


