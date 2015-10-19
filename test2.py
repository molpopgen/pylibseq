import pyseq.polytable as pypt


x = pypt.simData([ (0.1,"01"),(0.2,"10") ])

print x.size()

for i in range(2):
    print "element i is", x[i]
