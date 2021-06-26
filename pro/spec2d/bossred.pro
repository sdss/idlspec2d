pro bossred,planfile,flags,clobber

planfile=file_basename(planfile)
temp=strsplit(file_basename(planfile,'.par'),'-',/extract)
plate=temp[1]
mjd=temp[2]

cd,getenv('BOSS_SPECTRO_REDUX')+'/holtz/'+plate+'p'

;spreduce2d,planfile,/plates
;rm_combine_script,'spPlancomb-'+plate+'-'+mjd+'.par', /xyfit,/loaddesi,/plates, run2d='holtz'

spawn,['runrvs','spPlancomb-'+plate+'-'+mjd+'.par'],/noshell

end
