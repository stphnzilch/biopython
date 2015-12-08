import sys
import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC  #from biopython tutorial 3.3
from Bio import Restriction
from Bio.Restriction import Restriction_Dictionary


#name= raw_input("name sequance:  ") 
#raw_sequance = raw_input("enter sequance: ")
#raw_res = raw_input("enter restriction site(HindIII, EcoRI): ") 

import sys
import random
import argparse


def linebreak():
    print "\n"
    print "----------------------------"
    print "\n"
def GCcontent(G,C,frwd):
      G_C = 100 * float(G + C) / len(frwd)  
      print "GC content", G_C, "%"

def melttemp(G,C,A,T):
      MT = 64.9 + 41 *float( (G + C - 16.4) / (A + T + G + C))
      print  "melting temp", MT    
def main():
    seq = ''

    parser = argparse.ArgumentParser()
    parser.add_argument("bases")
    args = parser.parse_args()
    with open(args.bases, 'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
               
                name = (line.strip())
            else:
                seq += line.strip()

                
    just_restriction = Restriction.HindIII
    restriction_enzyme = just_restriction

    res_site =  Seq(restriction_enzyme.site)


    forward_resriction = res_site
    backward_restriction = res_site.reverse_complement()
    print just_restriction

    sequence_sites = forward_resriction + seq + backward_restriction

    print "restriction site", restriction_enzyme.site   #just the name of the input sequance

    linebreak()

    print name
    print "Length", len(sequence_sites)                 #  this section  shows the length and name of the input sequance
    linebreak()

    frwd = sequence_sites[0:20]   

    print  name, ", forward primer"  #prints the forward primer and the meling point


    print "forward", frwd

    G = frwd.count("G")
    C = frwd.count("C")
    A = frwd.count("A")
    T = frwd.count("T")

    GCcontent(G,C,frwd)

    melttemp(G,C,A,T)



    linebreak()       #prints the sequnce an meling point and gc content


    revs = sequence_sites[-20:]                                       #sequace for the reverse primer
    revs_revs_comp = revs.reverse_complement()                            #reverse complement of the last 20 characters
                                                                       #documentation same as forward

    print   name, ", reverse primer"
    print revs_revs_comp

    G = revs_revs_comp.count("G") 
    C = revs_revs_comp.count("C")
    A = revs_revs_comp.count("A")
    T = revs_revs_comp.count("T")
    GCcontent(G,C,frwd)
    melttemp(G,C,A,T)

main()        #EVERTHING!!!!!!
