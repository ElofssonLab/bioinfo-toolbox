import multiprocessing as mp
import time


def worker(number):
    print (str(number)+' start')
    for n in range(0, number, 1):
        time.sleep(2)
        print ('wait for '+str(number))
    print (str(number)+' done!')
    q.put(number)

def output():
    while 1: print (str(q.get())+' is a good boi!')

if __name__ == '__main__':


    joblist = [2, 4, 6, 8]

    p = mp.Process(target=output, args=())  ##### single process
    pool = mp.Pool(processes=4)             ##### set a pool of processes, jobs in a list are handled between those processes
    q = mp.Queue()                          ##### queue

    p.start()

    pool.map(worker, joblist)
    pool.close()
    pool.join()
    
    p.terminate()
