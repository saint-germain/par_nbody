import sys
import os
import numpy as np
# Folder where orbital data will be stored. This should not be changed,
# better to change wdir - working directory
odir = sys.argv[1]
odir += '/'

folder_to_save = sys.argv[2] + '/'

try:
    os.mkdir(odir)
    # print('Created folder: ' + odir)
except:
    pass
    # print('')

allfiles = os.listdir(odir)
txtFiles = [] # Lista para agregar los txt

for i in allfiles:
    txtFiles.append(odir + i) # Sistemas en txt

txtFiles = np.sort(txtFiles)

# outfile2 = open('systemList.txt','w')
# to_save = txtFiles[732:854]

# for texto in to_save:
#     outfile2.write(texto[11:28] + '\n')

# outfile2.close()
# exit()
# print(txtFiles[839])
# print(txtFiles[915])
# exit()

if __name__ == '__main__':

    from multiprocessing import Pool
    # from nbody6_SJM3 import integrator
    from mainIntegrator import mainIntegrator
    import parmap
    import datetime
    start = datetime.datetime.now()

    try:
        parmap.map(mainIntegrator, txtFiles[0:192], folder_to_save, pm_chunksize=1, pm_processes=4, pm_pbar=True)
    except Exception as error:
        print(error)

    # with parmap.map(mainIntegrator, txtFiles[610:615], pm_chunksize=1, pm_processes=4, pm_p) as result1:
        # data_task1 = None
        # task1_working = True
        # while task1_working:
            # if task1_working and result1.ready():
                # print("Task 1 has finished!")
                # data_task1 = result1.get()
                # task1_working = False

    # pool = Pool(3) # Como argumento, numero de procesos
    # pool.map(mainIntegrator, txtFiles[610:], 3)
    
    # pool.close()
    # pool.join()

    end = datetime.datetime.now()

    print('Global Timer:', end - start)
