from subprocess import call
import string

def plot_lambdaJ():

    lambdas = [0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1]

    name = 'test'
    length = 143
    pdb_acc = '../evfold/workingcopy/pdbs/LRR_5.pdb'
    pfam_acc = 'LRR_5_test'

    for l in lambdas:
        call(['Rscript', 'plot_contact_map.R', pdb_acc, pfam_acc + '_' + str(l), str(length), '../plmDCA_symmetric_v2/output/' + name + '_' + str(l) + '.scores', '0'])

def main():

    #pfam_ids = ["PF00023","PF00057","PF00090","PF00096","PF00191","PF00560","PF00626","PF00681","PF00683","PF01239","PF01391","PF01468","PF01473","PF01754","PF01826","PF01846","PF02412","PF02986","PF03989","PF05001","PF05552","PF06565","PF06898","PF07554","PF07661","PF09479"]

    pfam_hit_dict = {}

    f =  open('Pfam-hit_PDB.txt', 'r')

    line = f.readline()  ##skip file header
    line = f.readline()

    while not line == '':
        line_arr = string.split(line,'\t')

        pfam_acc = line_arr[0]
        pfam_hit_dict[pfam_acc] = line_arr[1:len(line_arr)]
    
        line = f.readline()

    print pfam_hit_dict

    for pfam_acc, prop_list  in pfam_hit_dict.iteritems():
        pfam_acc = pfam_acc.split('.')[0]
        print pfam_acc
        pdb_acc = prop_list[0].strip()
        chain = prop_list[1].strip()
        pos_start = prop_list[2].strip()
        pos_end = prop_list[3].strip()
        length = prop_list[6].strip()
    
        if pfam_acc[0:2] != 'CL':
            if pdb_acc != 'none':
                call(['Rscript', 'plot_contact_map.R', pdb_acc, pfam_acc, length, 'plmDCA/output/known_struct/' + pfam_acc + '_full.txt_scores.txt', pos_start])


def calculate_offset(pfam_acc, pdb_acc):
    
    


    return offset

if __name__ == "__main__":
    #main()
    plot_lambdaJ()
