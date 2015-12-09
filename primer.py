import sys
import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC  #from biopython tutorial 3.3
from Bio import Restriction
from Bio.Restriction import Restriction_Dictionary
import argparse
import primerlibrary

#linebreak prints a an empty line and a dotted line
def linebreak():
    print "\n"
    print "----------------------------"
    print "\n"


#GCcontent is a simple fraction that divides the amont of gs and cs and divides them by the whole length of the primer
def GCcontent(G,C,frwd):
      G_C = 100 * float(G + C) / len(frwd)  
      print "GC content", G_C, "%"


#meltingtemp uses a function that does syuff takes in a few variables and prints the thing
def melttemp(G,C,A,T):
      MT = 64.9 + 41 *float( (G + C - 16.4) / (A + T + G + C))
      print  "melting temp", MT

    
#countbase sets the global variable of GCTA and count them in the seq
def countbase(seq):
      global G
      global C
      global A
      global T
      G = seq.count("G")
      C = seq.count("C")
      A = seq.count("A")
      T = seq.count("T")

#restriction toss requires no input variables if the 
#restriction site needs to be changed it must be done in this function.
#first it defines the global variables
#then it names a variable the restriction enzyme pulling in 
#formation from a biopython dictionary
#then the name is turned into a the sequance of bases and converted to a string
#then the name and the seqaunce is printed.
#the second part does the same thing but resse complements the sequance
def restriction_toss():    
    global forward_res_site
    global backward_restriction
    global reverse_restriction_enzyme
    global forward_restriction_enzyme

    forward_res_name = Restriction.HindIII
    forward_res_site =  Seq(forward_res_name.site)
    print forward_res_name
    print "restriction site", forward_res_name.site   
   
    reverse_res_name = Restriction.EcoRI
    reverse_res_site =  Seq( reverse_res_name.site)
    backward_restriction = reverse_res_site.reverse_complement()
    print reverse_res_name
    print "restriction site",  reverse_res_name.site   
#divides the data coming in
#if the line stars with a > then its a titile
#then everything else is a base sequance
#store variable name and sequ
def parsing():
    global sequ
    global name
    global parser
    sequ = ''
    parser = argparse.ArgumentParser()
    parser.add_argument("bases")
    args = parser.parse_args()
    with open(args.bases, 'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
               
                name = (line.strip())
            else:
                sequ += line.strip()

#adds the restriction sites to the original sequance
def seq_adding():
    global sequence_sites
    sequence_sites = forward_res_site + sequ + backward_restriction



def main():
    parsing()
    restriction_toss()    
    seq_adding()

    linebreak()

    print name
    print "Length", len(sequence_sites)

    linebreak()

    frwd = sequence_sites[0:20]   
    seq = frwd
    print  name, ", forward primer"  
    print "forward", seq
    countbase(seq)
    GCcontent(G,C,seq)
    melttemp(G,C,A,T)

    linebreak()      

    revs = sequence_sites[-20:] #sequace for the reverse primer
    seq = revs.reverse_complement()#reverse complement of the last 20 characters                                         
    print   name, ", reverse primer"
    print "reverse", seq
    countbase(seq)
    GCcontent(G,C,seq)
    melttemp(G,C,A,T)

main()        #EVERYTHING!!!!!!
