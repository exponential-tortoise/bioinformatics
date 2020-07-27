import numpy as np 

def entropy(p):

    # p is a list
    ans=0
    for i in p:
        if i!=0:
            ans-=(i*np.log10(i))
    return ans

print(entropy([0.5,0,0,0.5]))
print(entropy([0.25, 0.25, 0.25, 0.25]))
print(entropy([0, 0, 0, 1]))
print(entropy([0.25, 0, 0.5, 0.25]))