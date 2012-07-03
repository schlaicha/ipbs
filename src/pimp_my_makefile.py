import sys
solvers=[ "BCGS_SSORk", "BCGS_NOPREC" , "CG_SSORk" , "CG_NOPREC" , "CG_Jacobi" , "CG_AMG_SSOR" , "BCGS_AMG_SSOR" ]
grids=[ "UGGRID" , "ALUGRID_SIMPLEX" ]
sys.stdout.write("EXTRA_PROGRAMS= ")
for i in range(len(solvers)):
    for j in range(1,4):
      for k in range(len(grids)):
        for l in range(2,4):
          sys.stdout.write("ipbs_"+grids[k]+"_"+str(l)+"d_"+solvers[i]+"_P"+str(j)+ " ")

sys.stdout.write("\n\n\n")

for i in range(len(solvers)):
    for j in range(1,4):
      for k in range(len(grids)):
        for l in range(2,4):
          sys.stdout.write("""
ipbs_"""+grids[k]+"""_"""+str(l)+"""d_"""+solvers[i]+"""_P"""+str(j)+"""_CPPFLAGS=$(shared_CPPFLAGS)""" \
    + """ -D"""+grids[k]+ """ -DGRIDDIM="""+str(l)+""" -DWORLDDIM="""+str(l)+
    """ -DPDEGREE="""+str(j)+""" -DLINEARSOLVER="""+str(i+1)+"""
ipbs_"""+grids[k]+"""_"""+str(l)+"""d_"""+solvers[i]+"""_P"""+str(j)+"""_SOURCES =$(ipbs_SOURCES)
ipbs_"""+grids[k]+"""_"""+str(l)+"""d_"""+solvers[i]+"""_P"""+str(j)+"""_LDFLAGS =$(ipbs_LDFLAGS)
        """)



