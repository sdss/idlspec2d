import numpy as np
import os
import sys
from astropy.io import fits as pyf
import importlib.machinery
import types
import os.path as ptt
location = os.getenv('IDLSPEC2D_DIR')
loader = importlib.machinery.SourceFileLoader('yanny', location+'/bin/yanny.py')
yanny = types.ModuleType(loader.name)
loader.exec_module(yanny)

redux=os.getenv('BOSS_SPECTRO_REDUX')
run2d=os.getenv('RUN2D')
run1d=os.getenv('RUN1D')
file_c=redux+'/'+run2d+'/conflist.fits'
hdu_list = pyf.open(file_c)
table_hdu = hdu_list[1]
table_data = table_hdu.data
plates=table_data.field('PLATE')
plates=np.unique(plates)
basedir=redux+'/'+run2d+'/'
ct_2d=0
cr_2d=0
cf_2d=0
ct_com=0
cr_com=0
ce_com=0
cf_com=0
cr_1d=0
for plate in plates:
  print('PLATE:'+str(plate))
  dir_p=basedir+str(plate).replace(' ','')+'p'+'/'
  os.system('ls '+dir_p+'spPlan2d* > list')
  f=open('list','r')
  platemjds_2d=[]
  for line in f:
     ct_2d=ct_2d+1
     file_2d=line.replace('\n','').replace(dir_p,'')
     log_2d=file_2d.replace('spPlan2d','spDiag2d').replace('.par','.log')
     try:
       f2=open(dir_p+'/'+log_2d,'r')
       lines=f2.readlines()
       last=lines[-1:][0].replace('\n','')
       #print(last)
       f2.close()
       cf_2d=cf_2d+1
       if not 'Successful completion' in last:
         print(file_2d)
         platemjd2d=file_2d.replace('spPlan2d-','').replace('.par','')
         platemjds_2d.extend([platemjd2d])
       else:
         #plan = yanny.read_yanny(dir_p+'/'+file_2d)
         #expo = plan['SPEXP']['name']
         #for itt in range(0, len(expo)):
         #   print(expo[itt][0].replace(''
         #sys.exit()
         cr_2d=cr_2d+1
     except:
       print(file_2d)
       platemjd2d=file_2d.replace('spPlan2d-','').replace('.par','')
       platemjds_2d.extend([platemjd2d])
  f.close()
 # sys.exit()
  os.system('ls '+dir_p+'spPlancomb* > list2')
  f=open('list2','r')
  for line in f:
     ct_com=ct_com+1
     file_2d=line.replace('\n','').replace(dir_p,'')
     log_2d=file_2d.replace('spPlancomb','spDiagcomb').replace('.par','.log')
     try:
     #aa=0
     #if aa == 0:
       f2=open(dir_p+'/'+log_2d,'r')
       lines=f2.readlines()
       last=lines[-1:][0].replace('\n','')
       #print(last)
       #print(file_2d)
       f2.close()
       cf_com=cf_com+1
       if not 'Successful completion' in last:
         print(file_2d,"T")
       else:
         plan = yanny.read_yanny(dir_p+'/'+file_2d)
         #print(plan)
         rawmjds = plan['SPEXP']['mjd']
         rawmjds = np.unique(rawmjds)
         t1=0
         for mjdc in rawmjds:
           for mjd2d in platemjds_2d:
              if str(mjdc).replace(' ','') in mjd2d:
                 #print(file_2d)
                 t1=1
         expo_dc=[]
         for linef in lines:
           if 'RM_SPCOMBINE_V5: WARNING: Discarding science exposure #' in linef:
              data=linef.replace('\n','').replace('RM_SPCOMBINE_V5: WARNING: Discarding science exposure #','').replace('\t','').split(' ')
              #print(data)
              #sys.exit()
              data=list(filter(None,data))
              expo_dc.extend([data[0]])
         if len(expo_dc) > 0:
           tpt=0
         else:
           tpt=1
         #print(expo_dc,tpt)
         #sys.exit()
         #a
         expo = plan['SPEXP']['name']
         for itt in range(0, len(expo)):
           file1t=expo[itt][0].replace('spFrame','spCFrame')
           file2t=expo[itt][1].replace('spFrame','spCFrame')
           if ptt.exists(dir_p+'/'+file1t) == False:
              if tpt==1:
                 t1=1
              else:
                 t1=1
                 #print("A",t1)
                 for expo_d in expo_dc:
                    if expo_d in file1t:
                      t1=0   
                 #print("B",t1)
           if ptt.exists(dir_p+'/'+file2t) == False:
              if tpt==1:
                 t1=1
              else:
                 t1=1
                 for expo_d in expo_dc:
                    if expo_d in file2t:
                      t1=0
         if t1 == 1:
           print(file_2d,'O')
         if t1 == 0:       
           cr_com=cr_com+1
         ce_com=ce_com+1
     except:
       tt=0
       print(file_2d,'Y')
     dir_1d=dir_p+'/'+run1d+'/'
     zline_1d=file_2d.replace('spPlancomb','spZline').replace('.par','.fits')
     if ptt.exists(dir_1d+zline_1d):
       cr_1d=1+cr_1d
     else:
       print(zline_1d)
     #   f3=open('list3','r')
     #linest=f3.readlines()
     #cr_1d=len(linest)+cr_1d
  f.close()
  #sys.exit()
print("Total of extracted plates-mjds: "+str(ct_2d))
print("Total of reduced extracted plates-mjs: "+str(cr_2d)+"/"+str(cf_2d))
print("Total of combined plates-mjds: "+str(ct_com))
print("Total of combined extracted plates-mjs: "+str(cr_com)+"/"+str(ce_com)+"/"+str(cf_com))
print("Total of 1d plates-mjs: "+str(cr_1d))
