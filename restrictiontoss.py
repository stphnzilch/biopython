import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC  #from biopython tutorial 3.3
from Bio import Restriction
from Bio.Restriction import Restriction_Dictionary

def restriction_toss():    
    global forward_resriction
    global backward_restriction
    just_restriction = Restriction.HindIII
    restriction_enzyme = just_restriction
    res_site =  Seq(restriction_enzyme.site)
    forward_resriction = res_site
    backward_restriction = res_site.reverse_complement()
    print just_restriction
    print "restriction site", restriction_enzyme.site   #just the name of the input sequance
    
def main():
     restriction_toss()
     print forward_resriction
     print backward_restriction

main()
