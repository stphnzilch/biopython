#G=0
#C=0
#A=0
#T=0

def GCcontent(G,C,seq):
      G_C = 100 * float(G + C) / len(seq)  
      print "GC content", G_C, "%"

def melttemp(G,C,A,T):
      MT = 64.9 + 41 *float( (G + C - 16.4) / (A + T + G + C))
      print  "melting temp", MT  
def countbase(seq):
      global G
      G = seq.count("G")
      global C
      C = seq.count("C")
      global A
      A = seq.count("A")
      global T
      T = seq.count("T")
def main():
    sequence_sites = "AAAGGATCTCTCTCTGAGAGAGAACCCCCTCTCTCTTAGAGAGACATCAT"
    revs = "GGATCTCTCTCTGAGAGAGAA"
    frwd = sequence_sites[0:20]   
    
    seq = frwd
    countbase(seq)
    print "frwd", G, A, T, C, len(frwd)
    GCcontent(G,C,seq)
    melttemp(G,C,A,T)

    seq = revs 
    countbase(seq)
    print "revs", G, A, T, C, len(revs)
    GCcontent(G,C,seq)
    melttemp(G,C,A,T)
main()
