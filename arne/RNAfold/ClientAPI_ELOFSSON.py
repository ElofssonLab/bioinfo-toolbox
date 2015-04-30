import random, sys, string, requests, json, hmac, hashlib, datetime, math, time

#You'll need a publicConsumerToken and a privateConsumerToken to use the API
publicConsumerToken = 'bdzswm4lvh'
privateConsumerToken = 'jn8lbgwzpg280wk126793oui5h'
url = 'https://api.denovodna.com/'

params = {}

def sendPost(methodName, params = {}):
    params['publicConsumerToken'] = publicConsumerToken
    params['timestamp'] = datetime.datetime.utcnow().isoformat()
    message = publicConsumerToken + params['timestamp'] + methodName
    params['signature'] = hmac.new(privateConsumerToken, message, hashlib.sha1).hexdigest()
    headers = {'Content-type': u'application/json', 'Accept': u'application/json'}
    try_counter = 0
    maxNumberAttempts = 5
    while True:
        r = requests.post(url + methodName,verify=False, data=json.dumps(params), headers=headers)
        output = r.json()
        if output['status'] == 'CREATED':
            break
        elif output['status'] == 'UNAUTHORIZED':
            try_counter += 1
            time.sleep(.5)
        elif try_counter > maxNumberAttempts:
            break
    return output

def sendGet(methodName, params = {}):
    params['publicConsumerToken'] = publicConsumerToken
    params['timestamp'] = datetime.datetime.utcnow().isoformat()
    message = publicConsumerToken + params['timestamp'] + methodName
    params['signature'] = hmac.new(privateConsumerToken, message, hashlib.sha1).hexdigest()
    headers = {'Content-type': u'application/json', 'Accept': u'application/json'}
    r = requests.get(url + methodName,verify=False, params=params, headers=headers)
#    print r
    return r.json()

def sendDelete(methodName, params = {}):
    params['publicConsumerToken'] = publicConsumerToken
    params['timestamp'] = datetime.datetime.utcnow().isoformat()
    message = publicConsumerToken + params['timestamp'] + methodName
    params['signature'] = hmac.new(privateConsumerToken, message, hashlib.sha1).hexdigest()
    headers = {'Content-type': u'application/json', 'Accept': u'application/json'}
    r = requests.delete(url + methodName,verify=False, params=params, headers=headers)
    return r.json()

def ReverseRBS(mRNA, startRange = None, rRNASeq = "acctcctta", algorithmVersion = "v1.1"):
    algorithm = 'ReverseRBS'
    params = {}
    params['mRNA'] = mRNA
    if startRange is None:
        params['start_range'] = [0, len(mRNA)]
    else:
        params['start_range'] = startRange
    params['16S_rRNA'] = rRNASeq
    params['algorithm_version'] = algorithmVersion
    params['algorithm'] = algorithm
    response = sendPost(algorithm, params)
    return response

def ForwardRBS(CDS, TargetTranslationInitRate = 1000, rRNASeq = "acctcctta", preSeq = "", RBSInit = None, maximizeTIR = False, algorithmVersion = "v1.1"):
    algorithm = 'ForwardRBS'
    params = {}
    params['CDS'] = CDS
    params['translation_init_rate'] = TargetTranslationInitRate
    params["16S_rRNA"] = rRNASeq
    params['pre_seq'] = preSeq
    params['RBS_init'] = RBSInit
    params['maximize_tir'] = maximizeTIR
    params['algorithm'] = algorithm
    params['algorithm_version'] = algorithmVersion
    response = sendPost(algorithm, params)
    return response

def ForwardRBSwithConstraints(CDS, TargetTranslationInitRate = 1000, RBSConstraints = "", rRNASeq = "acctcctta", preSeq = "", RBSInit = None, maximizeTIR = False, algorithmVersion = "v1.1"):
    algorithm = 'ForwardRBS'
    params = {}
    params['CDS'] = CDS
    params['translation_init_rate'] = TargetTranslationInitRate
    params['RBS_Constraints'] = RBSConstraints
    params["16S_rRNA"] = rRNASeq
    params['pre_seq'] = preSeq
    params['RBS_init'] = RBSInit
    params['maximize_tir'] = maximizeTIR
    params['algorithm'] = algorithm
    params['algorithm_version'] = algorithmVersion
    response = sendPost(algorithm, params)
    return response

def RBSLibraryCalculator_SearchMode(CDS, minTIR = "1", maxTIR = "100000", numVariants = "16", rRNASeq = "acctcctta", preSeq = "", RBSInit = None, algorithmVersion = "v1.0"):
    algorithm = 'RBSLibraryCalculator_SearchMode'
    params = {}
    params['CDS'] = CDS
    params['min_translation_init_rate'] = minTIR
    params['max_translation_init_rate'] = maxTIR
    params['expression_level_width'] = math.log10(maxTIR/minTIR) / int(numVariants)
    params["16S_rRNA"] = rRNASeq
    params['pre_seq'] = preSeq
    params['InitialRBS'] = RBSInit
    params['algorithm'] = algorithm
    params['algorithm_version'] = algorithmVersion
    params['verbose'] = True
    response = sendPost(algorithm, params)
    return response

def RBSLibraryCalculator_GenomeSearchMode(CDS, minTIR = "1", maxTIR = "100000", numVariants = "16", maxConsecutiveDegeneracy = "6", rRNASeq = "acctcctta", preSeq = "", RBSInit = None, algorithmVersion = "v1.0"):
    algorithm = 'RBSLibraryCalculator_GenomeSearchMode'
    params = {}
    params['CDS'] = CDS
    params['min_translation_init_rate'] = minTIR
    params['max_translation_init_rate'] = maxTIR
    params['MaxConsecutiveDegeneracy'] = maxConsecutiveDegeneracy
    params['expression_level_width'] = math.log10(maxTIR/minTIR) / int(numVariants)
    params["16S_rRNA"] = rRNASeq
    params['pre_seq'] = preSeq
    params['InitialRBS'] = RBSInit
    params['algorithm'] = algorithm
    params['algorithm_version'] = algorithmVersion
    params['verbose'] = True
    response = sendPost(algorithm, params)
    return response

def RBSLibraryCalculator_EvaluateMode(CDS, dRBSSequence, rRNASeq = "acctcctta", preSeq = "", algorithmVersion = "v1.0"):
    algorithm = 'RBSLibraryCalculator_GenomeSearchMode'
    params = {}
    params['CDS'] = CDS
    params["dRBS"] = dRBSSequence
    params['pre_seq'] = preSeq
    params["16S_rRNA"] = rRNASeq
    params['algorithm'] = algorithm
    params['algorithm_version'] = algorithmVersion
    params['verbose'] = True
    response = sendPost(algorithm, params)
    return response

#Available at a later time.
def OperonCalculator_EvaluateMode(promoter = "", mRNA = "", rRNASeq = "acctcctta", algorithmVersion = "v1.0"):
    algorithm = 'OperonCalculator_EvaluateMode'
    params = {}
    params['promoter' ] =promoter
    params['mRNA'] = mRNA
    params['16S_rRNA'] = rRNASeq
    params['algorithm'] = algorithm
    params['algorithm_version'] = algorithmVersion
    response = sendPost(algorithm, params)
    return response

def getJobStatuses():
    jobList = sendGet('JobStatus', {})
    return jobList

if __name__ == "__main__":
    #Examples of different types of inputs
    mRNA = "".join([random.choice(('A','G','C','T')) for x in range(35)]) + 'ATG' + "".join([random.choice(('A','G','C','T')) for x in range(35)])
    CDS = 'ATG' + "".join([random.choice(('A','G','C','T')) for x in range(50)])
    targetTIR = 10**(random.random() * 5)
    RBSConstraints = 'TCTAGA' + 'N'*25 + 'GAATTC'
    minTIR = 10**(random.random() * 2)
    maxTIR = minTIR + 10**(3 + random.random() * 2)
    numVariants = 16
    dRBS = [random.choice(('A','G','C','T')) for x in range(35)]
    for i in range(3):
        pos = random.randint(0,34)
        dRBS[pos] = random.choice( ('W','S','M','K','R','Y','B','D','H','V','N') )
    dRBS = "".join(dRBS)
    promoter = "".join([random.choice(('A','G','C','T')) for x in range(20)])
    long_mRNA = "".join([random.choice(('A','G','C','T')) for x in range(100)]) + 'ATG' + "".join([random.choice(('A','G','C','T')) for x in range(999)]) + 'ATG' + "".join([random.choice(('A','G','C','T')) for x in range(9)]) + 'TAA' +  "".join([random.choice(('A','G','C','T')) for x in range(999)]) + 'TAA'
    sequence=[]
    list=[]
    expression={}
    with open(sys.argv[1]) as fp:
        for line in fp:
            sequence=line.strip().split(' ')
            list.append(sequence[1])
            expression[sequence[1][0:54]]=sequence[2]
    fp.close()
#How to Call Jobs
    for seq in list:
        mRNA=seq
        response = ReverseRBS(mRNA=mRNA)
#        print response
        #How to Retrieve Results
    resultList= sendGet('Result', {})
    reply={}
    for result in resultList:
        status='TEST'
        while status != 'Finished': 
            response = sendGet('Result',{'id' : result.split('/')[-1]})
            for key,value in response.iteritems():
#                print key,value
                reply[key]=value
                if key=='status':
                    status=value
                if type(value)==type(dict()):
                    for keyB,valueB in value.iteritems():
                        reply[key][keyB]=valueB
                        if keyB=='status':
                            status=valueB
                        if type(valueB)==type(dict()):
                            for keyC,valueC in value.iteritems():
                                reply[key][keyB][keyC]=valueC
                                if keyC=='status':
                                    status=valueC
                        if valueB=='RBS_list':
                            for keyC,valueC in value.iteritems():
                                print "test: ",keyC,valueC

        j=0
        for i in reply['ReverseRBS']['RBS_list']:
            if reply['ReverseRBS']['RBS_list'][j]['start_position']==20:
                print reply['ReverseRBS']['RBS_list'][j]['sequence'],reply['ReverseRBS']['RBS_list'][j]['start_position'],reply['ReverseRBS']['RBS_list'][j]['tir'],expression[reply['ReverseRBS']['RBS_list'][j]['sequence'][0:54]]
#            print 'mRNA:  ',reply['ReverseRBS']['RBS_list'][j]['sequence']
#            print 'Start: ',reply['ReverseRBS']['RBS_list'][j]['start_position']
#            print 'TIR: ',reply['ReverseRBS']['RBS_list'][j]['tir']
            j=j+1
#        print "sendGET: ",json.dumps(response,indent=1)
        response = sendDelete('Result', {'id' : result.split('/')[-1]})
#        print "sendDelete: ",response


#AraH: TAATTTTGTTTAACTTTAAGAAGGAGACNNNNNNATGNNNNNNTCTGTTTCTACATCGGGGTCTGGCGCACCTAAGTCGTCATTCAGCTTCGGGCGTATCTGGGATCAGTACGGCATGCTGGTGGTGTTTGCGGTGCTCTTTATCGCCTGTGCCATTTTTGTCCCAAATTTTGCCACCTTCATT     (ATG pos=31)
#                 TTTAAGAAGGAGACCTGAGTATGATGCGTTCTGTTTCTACATCGGGGTCTGGCGCAC  (ATG pos=20)

#NarK: TAATTTTGTTTAACTTTAAGAAGGAGACNNNNNNATGNNNNNNTCATCCGCCCCCGAAAGGGCTACTGGAGCTGTCATTACAGATTGGCGACCGGAAGATCCTGCGTTTTGGCAACAACGCGGTCAACGTATTGCCAGCCGCAACCTGTGGATTTCCGTTCCCTGTCTGCTGCTGGCGTTTTGC     (ATG pos=31)
#                 TTTAAGAAGGAGACTTGTATATGTCGCATCACTCATCCGCCCCCGAAAGGGCTACTG  (ATG pos=20)
