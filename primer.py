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
print "  "
print "--------------------------------------------------"
print "  "

print name
print "Length", len(sequence_sites)                 # just the input sequancse                 
print "  "
print "--------------------------------------------------"
print "  "
frwd = sequence_sites[0:20]   

print  name, ", forward primer"  #frwd_raw_res       print ">forward", name , #frwd_raw_


print frwd

frwd_G = frwd.count("G")
frwd_C = frwd.count("C")
frwd_A = frwd.count("A")
frwd_T = frwd.count("T")


frwdG_C = 100 * float(frwd_G + frwd_C) / len(frwd)  

print "GC content", frwdG_C, "%"

frwd_MT = 64.9 + 41 *float( (frwd_G + frwd_C - 16.4) / (frwd_A + frwd_T + frwd_G + frwd_C))

print  "melting temp", frwd_MT  

print "  "
print "----------------------------------------------------"
   

print"  "


revs = sequence_sites[-20:]                                       #sequace for the reverse primer
revs_revs_comp = revs.reverse_complement()                            #reverse complement of the last 20 characters
                                                                       #documentation same as forward

print   name, ", reverse primer"
print revs_revs_comp

revs_G = revs_revs_comp.count("G") 
revs_C = revs_revs_comp.count("C")
revs_A = revs_revs_comp.count("A")
revs_T = revs_revs_comp.count("T")
revsG_C = 100 * float (revs_G + revs_C ) / len(revs_revs_comp)
print "GC content: ",revsG_C, "%"
revs_MT =  64.9 + 41 *float( (revs_G + revs_C - 16.4) / (revs_A + revs_T + revs_G + revs_C))    #http://www.biophp.org/minitools/melting_temperature/demo.php?formula=basic
print "melting temp",  revs_MT
