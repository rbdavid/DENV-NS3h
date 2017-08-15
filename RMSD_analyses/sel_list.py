#!/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python
##!/Library/Frameworks/Python.framework/Versions/2.7/bin/python

sel = []
### Paper RMSD selections ### 
sel.append(['a2_subdomain1_backbone','backbone and resid 57:68 and not name H*'])
sel.append(['motif_2_backbone','backbone and resid 117:124 and not name H*'])
sel.append(['aligned_CAs','protein and (resid 20:25 50:55 73:75 90:94 112:116 142:147 165:169 190:194 214:218 236:240 253:258 303:307) and name CA'])
sel.append(['aligned_betas','protein and (resid 20:25 50:55 73:75 90:94 112:116 142:147 165:169 190:194 214:218 236:240 253:258 303:307) and not name H*'])

#sel.append(['',''])
