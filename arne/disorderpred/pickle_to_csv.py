#/usr/bin/env python3



with open('outpred_rna.pickle','rb') as f:
    rnadata = pickle.load(f)

with open('outpred_pro.pickle','rb') as f:
    prodata = pickle.load(f)

with open('mobidata_K.pickle','rb') as f:
    mobidata = pickle.load(f)



with open(self.file_path, 'a') as csv_file:
        wr = csv.writer(csv_file, delimiter='|')
        row = ['klines', symbol]
        row.extend(pickle.dumps(object))
        wr.writerow(row)

    
