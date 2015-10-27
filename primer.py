
import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC  #from biopython tutorial 3.3
from Bio import Restriction
from Bio.Restriction import Restriction_Dictionary

name= raw_input("name sequance:  ") 

raw_sequance = raw_input("enter sequance: ") 
  

rwd_raw_res = raw_input("enter restriction site: ") 

print rwd_raw_res 

rest_enzyme = Restriction.HindIII

sequance_sites = raw_sequance

print "restriction site", rest_enzyme.site

print "Length", len(sequance_sites)                                   

frwd = sequance_sites[0:20]   

print ">forward", name , #frwd_raw_res       print ">forward", name , #frwd_raw_r




print frwd

frwd_G = frwd.count("G")
frwd_C = frwd.count("C")
frwd_A = frwd.count("A")
frwd_T = frwd.count("T")


frwdG_C = 100 * float(frwd_G + frwd_C) / len(frwd)  

print "GC content", frwdG_C, "%"

frwd_MT = 64.9 + 41 *float( (frwd_G + frwd_C - 16.4) / (frwd_A + frwd_T + frwd_G + frwd_C))

print  "melting temp", frwd_MT  



revs = Seq(sequance_sites[-20:])                                       #sequace for the reverse primer
revs_revs_comp = revs.reverse_complement()                            #reverse complement of the last 20 characters
                                                                       #documentation same as forward

print  ">reverse", name 
print revs_revs_comp

revs_G = revs_revs_comp.count("G") 
revs_C = revs_revs_comp.count("C")
revs_A = revs_revs_comp.count("A")
revs_T = revs_revs_comp.count("T")
revsG_C = 100 * float (revs_G + revs_C ) / len(revs_revs_comp)
print "GC content: ",revsG_C, "%"
revs_MT =  64.9 + 41 *float( (revs_G + revs_C - 16.4) / (revs_A + revs_T + revs_G + revs_C))    #http://www.biophp.org/minitools/melting_temperature/demo.php?formula=basic
print "melting temp",  revs_MT
