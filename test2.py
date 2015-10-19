import libsequence.polytable as pypt


x = pypt.simData()
x.assign([ (0.1,"01"),(0.2,"10") ])

print x.size()

for i in range(2):
    print "element i is", x[i]

pos = [0.1,0.2,0.3,0.4]
data = ["0101","1011"]

x.assign_sep(pos,data)

for i in range(2):
    print "element i is", x[i]

#x = pypt.polyTable()

print x.size()
