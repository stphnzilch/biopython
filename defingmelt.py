
def GCcontent(G,C,frwd):
      G_C = 100 * float(G + C) / len(frwd)  
      print "GC content", G_C, "%"

def melttemp(G,C,A,T):
      MT = 64.9 + 41 *float( (G + C - 16.4) / (A + T + G + C))
      print  "melting temp", MT  
def main():
    sequence_sites = "AAAGGATCTCTCTCTGAGAGAGAACCCCCTCTCTCTTAGAGAGACATCAT"

    frwd = sequence_sites[0:20]   
    G_C = 0
    MT = 0

    G = frwd.count("G")
    C = frwd.count("C")
    A = frwd.count("A")
    T = frwd.count("T")

    GCcontent(G,C,frwd)

    melttemp(G,C,A,T)
    print G, A, T, C, len(frwd)
   
    

main()
