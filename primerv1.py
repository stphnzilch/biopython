


from Bio.Seq import Seq
from Bio.Alphabet import IUPAC  #from biopython tutorial 3.3
from Bio import Restriction
from Bio.Restriction import Restriction_Dictionary

name= raw_input("name sequance:  ")                                    #input a name
raw_sequance = raw_input("enter sequance: ")                           #input a sequance you want a primer for
#frwd_raw_res = HindIII      #raw_input("enter the forward restriction site: ")       #input the priimer you want on the forward sidee
#revs_raw_res = EcoRI      # raw_input("enter the reverse restriction site: ")       #input the primer you want on the reverse side

#from Bio.Restriction import (frwd_raw_res)
res = Restriction.EcoRI.site
#res = Restriction(frwd_raw_res)
print res




#print  Restriction.HindIII.site                                  #!!!!!!problem prints the restricion site eg. for EcoRI GAATTC
#print  Restriction.EcoRI.site                                   #!!!!!!problem prints the restricion site eg. for EcoRI GAATTC
print ">", name,# HindIII, EcoRI                            # prints the Fasta tag name + forward resticion enzyme + reverse restricion enzyme
#sequance_sites =  frwd_raw_res.site  + raw_sequance + revs_raw_res.site #!!!!!problem makes a variable of the full
                                                                         #sequance qith the forward and reverse restricion sites 

sequance_sites = raw_sequance                                          #while the previous line is broken this makes the program
                                                                       # work by simply printing  the raw sequance without the restricion sites
print "Length", len(sequance_sites)                                   #prints the working sequance (eventully with the res sites)

frwd = sequance_sites[0:20]                                           #forward primer-makes a variable of the first 20 charictars 
#print "length forward primer", len(frwd)                             #prints the length of the forward pirmer 20 duh

print ">forward", name , #frwd_raw_res                                    #prints the fasta tag for forwward primer with name of sequance and the cuting site
print frwd                                                            #prints the forward primer (eventully with the primer)
frwd_G = frwd.count("G")                                              #counts G's
#print "G: ", frwd_G                                                  #prints number of G's
frwd_C = frwd.count("C")                                              #counts the number of C's
#print "C: ", frwd_C                                                  #prints number of C's
frwd_A = frwd.count("A")                                              #counts number of A's
frwd_T = frwd.count("T")                                              #counts the number of T's
frwdG_C = 100 * float(frwd_G + frwd_C) / len(frwd)                    # GC content (G/C)/len*100
print "GC content", frwdG_C, "%"                                      #prints the gc contents in a precent
frwd_MT = 64.9 + 41 *float( (frwd_G + frwd_C - 16.4) / (frwd_A + frwd_T + frwd_G + frwd_C))   #http://www.biophp.org/minitools/melting_temperature/demo.php?formula=basic
                                                                      #the melting temp equation
print  "melting temp", frwd_MT                                        #prints the melting temp

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
