import numpy as np
from scipy.special import comb
import itertools
import sys
n=int(sys.argv[1])
def buildbin(n):
    if n==0:
        return np.zeros((0,1),dtype=np.uint8)
    else:
        m=buildbin(n-1)
        M=np.zeros((n,2**n),dtype=np.uint8)
        M[0:(n-1),0:2**(n-1)]=m
        M[0:(n-1),2**(n-1):2**n]=m
        M[n-1,0:2**(n-1)]=0
        M[n-1,2**(n-1):2**n]=1
        return M

if __name__ == '__main__':
    f=open('tetra_forbidden.lp','w')
    var_num=2**(n-1)-1+4*comb(n,4,exact=True)
    constraint_num=4*comb(n,4,exact=True)*(2**(n-4))+comb(n,4,exact=True)
    var_name_arr=['x'+str(i) for i in range(var_num)]
    f.write('Maximize\n')
    f.write('obj: '+' + '.join(var_name_arr[4*comb(n,4,exact=True):])+'\n')
    f.write('Subject To\n')

    #produce constraints
    #first 4C(n,4) variebles are local splits
    local_id=np.asarray(list(itertools.combinations(range(n),4)))
    local_offset=4*comb(n,4,exact=True)
    constraint_id=0
    for i in range(comb(n,4,exact=True)):
        constraint_id+=1
        f.write(f' c{constraint_id}: ')
        f.write(' + '.join(['x'+str(4*i+k) for k in range(4)]))
        f.write(' <= 3\n')
    #relax var constraints
    bin_ind=np.transpose(buildbin(n-4))
    for i in range(comb(n,4,exact=True)):
        s=local_id[i,:]
        res=list(set(range(n))-set(s))
        for j in range(bin_ind.shape[0]):
            for k in range(4):
                split_list=[s[k]]+[res[i] for i in range(n-4) if bin_ind[j,i]==1]
                #print(split_list)
                split_ind=np.sum(2**np.array(split_list))
                if n-1 in split_list:
                    split_ind=2**n-1-split_ind
                split_ind-=1
                split_ind+=local_offset
                constraint_id+=1
                f.write(f' c{constraint_id}: x{split_ind} - x{4*i+k} <= 0\n')
    f.write('Bounds\nBinary\n')
    f.write(' '.join(var_name_arr)+'\n')
    f.write('End\n')
    print(var_num,constraint_num)
    f.close()
