
def parse(afile, sep=','):
    print "sep=" + sep
    contacts = []

    for aline in afile:
        if aline.strip() != '':
            line_arr = aline.strip().split(sep)
            if line_arr[0].startswith('E'):
                continue
            i = int(line_arr[0])
            j = int(line_arr[1])
            score = float(line_arr[-1])
            if abs(i - j) > 4:
                contacts.append((score, i, j))

    contacts.sort(key=lambda x: x[0], reverse=True)
    #print contacts
    return contacts
    #return [contacts_x, contacts_y, scores] 


def write(contacts, afile, sep=','):
    
    for c in contacts:
        line = '%s%s%s%s%s\n' % (c[1], sep, c[2], sep, c[0])
        afile.write(line)
