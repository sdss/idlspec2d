;+
; NAME:
;   match_tsobjtargets
; PURPOSE:
;   Match a set of ra, dec objects to the tsobjtargets files
; COMMENTS:
; CALLING SEQUENCE:
;   match_tsobjtargets,tsobjtargetsbase,ra,dec,run,rerun,camcol,field,id
; INPUTS:
;   tsobjtargetsbase - base directory for tsobjtargets to match to
;   ra,dec - equatorial coords (in degrees)
;   run, rerun, camcol, field, id - SDSS id number
; OPTIONAL INPUTS:
;   matchlength - in arcsec, how close a match needs to be [1.]
;   tsobjbase - for objects which are not in tsobjtargets, this is the
;               base directory
; KEYWORDS:
; OUTPUTS:
;   tsobj - tsobj entries found (all zeroes if none)
; OPTIONAL OUTPUTS:
; PROCEDURES CALLED:
; BUGS:
;   Fails for rerun numbers > 1000
; REVISION HISTORY:
;   2002-Jan-11  Written by Blanton (NYU)
;-
pro match_tsobjtargets,tsobjtargetsbase,ra,dec,run,rerun,camcol,field,id,tsobj,matchlength=matchlength,tsobjbase=tsobjbase

maxrerun=1000l
maxfield=2000l

if (NOT keyword_set(matchlength)) then matchlength=1.d
matchlength=matchlength/3600.

tmpindx=run*maxrerun*6l+rerun*6l+camcol-1l
sortindx=tmpindx[sort(tmpindx)]
indx=long(sortindx[uniq(sortindx)])

for i=0l, n_elements(indx)-1l do begin
; get the file name, if it exists
    curr_run=indx[i]/(maxrerun*6l)
    curr_rerun=(indx[i]-curr_run*maxrerun*6l)/6l
    curr_camcol=indx[i]-curr_run*maxrerun*6l-curr_rerun*6l+1l
    curr_run5=string(curr_run,format='(i5.5)')
    curr_rerun0=string(curr_rerun,format='(i0)')
    curr_camcol1=string(curr_camcol,format='(i1)')
    searchpath=tsobjtargetsbase+'/runs/'+curr_run5+'/'+ $
                      curr_rerun0+'/tsObjTargets-'+curr_run5+'-'+ $
                      curr_rerun0+'-'+curr_camcol1+'-*'
    splog,'searching '+searchpath
    filename=findfile(searchpath)
    if (filename[0] ne '') then begin
        for j=0l, n_elements(filename)-1l do begin
            splog,filename[j]
; find tsobj entries for each
            tmptsobj=mrdfits(filename[j],1)
            if (keyword_set(tmptsobj)) then begin
                if (n_elements(tmptsobj) gt 0) then begin
                    match1=lonarr(1)-1l
                    spherematch,ra,dec,tmptsobj.ra,tmptsobj.dec,matchlength, $
                      match1,match2,distance12
                    if (match1[0] ne -1) then begin
                        if (NOT keyword_set(tsobj)) then begin
                            dum={dummy, blah:0.}
                            outblank=tmptsobj[match2[0]]
                            struct_assign,dum,outblank
                            tsobj=replicate(outblank,n_elements(ra))
                        endif
                        tsobj[match1]=tmptsobj[match2]
                    endif
                endif
            endif
        endfor
    endif
endfor

; Check if anything is unmatched
moretodo=0l
if(NOT keyword_set(tsobj)) then begin
    blankindx=lindgen(n_elements(ra))
    moretodo=1l
endif else begin
    blankindx=where(tsobj.run eq 0,count)
    if(count gt 0) then moretodo=1l
endelse

; If they are, and we know where to look for the full tsobj files,
; try to find them there
if(moretodo gt 0l and keyword_set(tsobjbase)) then begin
    tmpindx=run[blankindx]*maxrerun*6l*maxfield+rerun[blankindx]*6l*maxfield $
      +(camcol[blankindx]-1l)*maxfield+field[blankindx]
    sortindx=tmpindx[sort(tmpindx)]
    indx=long(sortindx[uniq(sortindx)])
    for i=0l, n_elements(indx)-1l do begin
; get the file name, if it exists
        curr_run=indx[i]/(maxrerun*6l*maxfield)
        curr_rerun=(indx[i]-curr_run*maxrerun*6l*maxfield)/(6l*maxfield)
        curr_camcol=(indx[i]-curr_run*maxrerun*6l*maxfield $
                     -curr_rerun*6l*maxfield)/maxfield+1l
        curr_field=(indx[i]-curr_run*maxrerun*6l*maxfield $
                    -curr_rerun*6l*maxfield-(curr_camcol-1l)*maxfield)
        curr_run6=string(curr_run,format='(i6.6)')
        curr_rerun0=string(curr_rerun,format='(i0)')
        curr_camcol1=string(curr_camcol,format='(i1)')
        curr_field4=string(curr_field,format='(i4.4)')
        searchpath=tsobjbase+'/'+curr_run6+'/'+ $
          curr_rerun0+'/tsObj-'+curr_run6+'-'+ $
          curr_camcol1+'-'+curr_rerun0+'-'+curr_field4+'.fit*'
        splog,'searching '+searchpath
        filename=findfile(searchpath)
        if (filename[0] ne '') then begin
            for j=0l, n_elements(filename)-1l do begin
                splog,filename[j]
; find tsobj entries for each
                tmptsobj=mrdfits(filename[j],1)
                if (keyword_set(tmptsobj)) then begin
                    if (n_elements(tmptsobj) gt 0) then begin
                        match1=lonarr(1)-1l
                        spherematch,ra[blankindx],dec[blankindx],tmptsobj.ra, $
                          tmptsobj.dec,matchlength,match1,match2,distance12
                        if (match1[0] ne -1) then begin
                            if (NOT keyword_set(tsobj)) then begin
                                dum={dummy, blah:0.}
                                outblank=tmptsobj[match2[0]]
                                struct_assign,dum,outblank
                                tsobj=replicate(outblank,n_elements(ra))
                            endif
                            tsobj[blankindx[match1]]=tmptsobj[match2]
                        endif
                    endif
                endif
            endfor
        endif
    endfor
endif

end

