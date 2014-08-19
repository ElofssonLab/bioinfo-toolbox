import sys
sys.path.append("/home/mircomic/bioinfo-toolbox")
from parsing import parse_contacts

# command line input
infile_name = sys.argv[1]

# guessing separator of constraint file
test_line = open(infile_name,'r').readline()
if len(test_line.split(',')) != 1:
    sep = ','
elif len(test_line.split(' ')) != 1:
    sep = ' '
else:
    sep = '\t'

# parse constraint file
c_list = parse_contacts.parse(open(infile_name, 'r'))
    
# sort contacts and write simple string
for c in c_list:
    print ('%s %s %s' % (c[1], c[2], c[0]))
