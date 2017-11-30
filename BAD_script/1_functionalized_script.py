#This Crosslink script is for updating the bond, angle, and dihedral from initial bond
#to obtain a functionalized COOH + gnp.
#Created 2 Mar 2017 (Choi)
#MEEM 1400 Michigan Technological University
#Houghton, MI, 49931

#Edited finding all angle & dihedral connections and all their types, 10 Mar 2017 (Choi)
#Tested with EEC molecule to check number of angle and dihedral that were generated and 
#to confirm epoxdide group, 12 Mar 2017
#Edited assigning angle_coeffs and dihedral_coeffs comment, 14 Mar 2017 (Choi)
#Edited add new angle and dihedral types based on BAD and vdW no., 23 Mar 2017 (Choi)
#Edited read and assign oplsaa forcefield file to assign angle and dihedral coefficients ...
#and get inputfile and name outputfile, 03 Apr 2017 (Choi)
#Edited fixing determining dihedral function to get all dihedral, 09 Apr 2017 (Choi)
#Edited fixing determining dihedral function back and forth, 25 May 2017 (Choi) line 270 delete"-1"

print '\n'
print 'Functionalized update starts now ~~\n'

import read1
import write
import sys
import operator
import time

class Bond: pass #Class for creating bond's objects
class Angle: pass #Class for creating angle's objects
class Dihedral: pass #Class for creating dihedral's object
class Ref: pass #Class for creaing reference's object
class Libatom: pass #Class for library of atoms' connection
class Libangle: pass #Class for library of angles' types
class Libdih: pass #Class for library of dihedral's types

filename = sys.argv[1]
m = read1.Molecule_File(filename)
ff = read1.Forcefield("oplsaa-tinker.prm")

#assign all oplsaa force field parameters into the variables
aff = ff.angleff
dff = ff.dihedralff
naff = ff.naff
ndff = ff.ndff

a = m.atoms
bnew = m.bonds
nnew = m.angles
dnew = m.dihedrals

b = bnew.copy() #reserved the new dictionary for updating all changed bonds
n = nnew.copy() #reserved the new dictionary for updating all changed angles
d = dnew.copy() #reserved the new dictionary for updating all changed dihedrals

mass = m.masses

natom = m.natomtypes
nbond = m.nbonds
nangle = m.nangles
ndihe = m.ndihedrals

m.extrabondperatom = 2 #extra bond per atom is necessary for new version of lammps

#--------------------------------------------------------------------------'
#----------------------vdW-&-BAD-number-in-OPLSAA--------------------------'
#--------------------------------------------------------------------------'

vdw = {} # list vdw number in OPLSAA
BAD = {} # list BAD (Bond-Angle-Dihedral) number in OPLSAA 
BAD[0] = 0
vdw[0] = 0
i = 1
pc = m.pair_coeffs
while i <= natom:
 nvdw = pc[i][2].split("BAD")[0].split("vdw")[1]
 nbad = pc[i][2].split("BAD")[1]
 vdw[i] = int(nvdw)
 BAD[i] = int(nbad)
 print vdw[i], BAD[i]
 i= i+ 1


#--------------------------------------------------------------------------'
#--------------------------Crosslink-information---------------------------'
#--------------------------------------------------------------------------'

#bini = 5 #initail bond type from  bond/create
#Cgnp_old = 1 #first initial old carbon atom type form gnp
#Cgnp_new = 5 #first initial new carbon atom type form gnp
#Ccarb_old = 2 #first initial old carbon atom type form carboxyl
#Ccarb_new = 6 #first initial new carbon atom type form carboxyl

#bgnp = 1 #bond type in pristine gnp [sp2-sp2]
#bsp23 = 6 #bond type between carbon sp2 and carbon sp3 [sp2-sp3]
#bcc_o = 4 #bond (double) C=O in carboxyl group 
#bcco = 2 #bond C-O in carboxyl group
#bcoh = 3 #bond O-H in carboxyl group

#nccc = 3 #angle type of Cgnp[old] - Cgnp[new] - Ccarb

#--------------------------------------------------------------------------'
#-----------------------------Confirm-atom-type----------------------------'
#--------------------------------------------------------------------------'

# Confirm atom types in the reactions
#for bidc in b:
# if b[bidc].type == bini:
#  if a[b[bidc].atomids[0]].type == Ccarb_old:
#        a[b[bidc].atomids[0]].type = Ccarb_new
#  elif a[b[bidc].atomids[1]].type == Ccarb_old:
#        a[b[bidc].atomids[1]].type = Ccarb_new
#  elif a[b[bidc].atomids[0]].type == Cgnp_old:
#        a[b[bidc].atomids[0]].type = Cgnp_new
#  elif a[b[bidc].atomids[1]].type == Cgnp_old:
#        a[b[bidc].atomids[1]].type = Cgnp_new
#  print a[b[bidc].atomids[0]].type, a[b[bidc].atomids[1]].type, '\n'

#--------------------------------------------------------------------------'
#--------------------------------Change-bond-------------------------------'
#--------------------------------------------------------------------------'

# Change bond to another bond
#print "____________________________________________________________________\n"
#print "Change to new bonds\n"

#for bid0 in b:
# if b[bid0].type == bini:
#  if a[b[bid0].atomids[0]].type == Ccarb_new:
#   Cgnp = b[bid0].atomids[1]
#   Ccarb = b[bid0].atomids[0]
#  elif a[b[bid0].atomids[0]].type == Cgnp_new:
#   Cgnp = b[bid0].atomids[0]
#   Ccarb = b[bid0].atomids[1]
#  print Cgnp , Ccarb
#  for bid in b:
#    if b[bid].type == bgnp:
#      if (a[b[bid].atomids[0]].type == Cgnp_new) or (a[b[bid].atomids[1]].type == Cgnp_new):
#        bnew[bid].type = bsp23
#        print 'Bond ID:', bid, 'was changed to bond type ==> ', bnew[bid].type
#        print b[bid].atomids[0], b[bid].atomids[1]
#    continue

#--------------------------------------------------------------------------'
#------------------------------Create-library--------------------------------'
#--------------------------------------------------------------------------'
print "____________________________________________________________________\n"
print "Create library\n"
liball = {}
for aid in a:
  countll = 0
  ll = {}
  ll[0] = 9999999
  ll[1] = 9999999
  ll[2] = 9999999
  ll[3] = 9999999
  for bid in b:
    if b[bid].atomids[0] == aid:
      ll[countll] = b[bid].atomids[1]
     # print aid, ll[countll]
      countll = countll + 1
    elif b[bid].atomids[1] == aid:
      ll[countll] = b[bid].atomids[0]
     # print aid, ll[countll]
      countll = countll + 1
 # print ll
  s_ll = sorted(ll, key=ll.get)
 # print  ll[s_ll[0]], ll[s_ll[1]], ll[s_ll[2]], ll[s_ll[3]]
  if (ll[s_ll[1]] == 9999999) and (ll[s_ll[2]] == 9999999) and (ll[s_ll[3]] == 9999999):
    sorted_atoms = [ll[s_ll[0]]]
  elif (ll[s_ll[1]] != 9999999) and (ll[s_ll[2]] == 9999999) and (ll[s_ll[3]] == 9999999):
    sorted_atoms = [ll[s_ll[0]],ll[s_ll[1]]]
  elif (ll[s_ll[1]] != 9999999) and (ll[s_ll[2]] != 9999999) and (ll[s_ll[3]] == 9999999):
    sorted_atoms = [ll[s_ll[0]],ll[s_ll[1]],ll[s_ll[2]]]
  elif (ll[s_ll[1]] != 9999999) and (ll[s_ll[2]] != 9999999) and (ll[s_ll[3]] != 9999999):
    sorted_atoms = [ll[s_ll[0]],ll[s_ll[1]],ll[s_ll[2]],ll[s_ll[3]]]
  #print "1", sorted_atoms  
  lll = Libatom()
  lll.atomids = sorted_atoms  
  liball[aid] = lll
  print aid, liball[aid].atomids
#sys.exit()
##--------------------------------------------------------------------------'
##------------------------------Create-angle--------------------------------'
##--------------------------------------------------------------------------'
#print "____________________________________________________________________\n"
#print "Create new angles\n"
b = bnew.copy()
l = liball # atom connection library
nnnew = {} #new library for all angles
libangle = {} # angle type library for collecting all angle types
nnn = Libangle()
nnn.atomtypes = [0,0,0] # angle type zero using as the dummy for first checking
libangle[0] = nnn
alibtype = 1 #count or id of angle type
count0 = 1
ncount = 1 #count number of angle in datafile
for lid in l:
  numl = len(l[lid].atomids) 
  count1 = 0
  while count1 < numl - 1:
    nl = l[lid].atomids[count1]
    nc = lid
    count2 = count1 + 1
    while count2 < numl:
      nr = l[lid].atomids[count2]
      nlt = a[nl].type
      nct = a[nc].type
      nrt = a[nr].type
      #print count0, nl, nc, nr
      #print 'type', nlt, nct, nrt
      count0 = count0 + 1
      count2 = count2 + 1
      # Check whether the angle type has already existed or not?
      check = 0 # checking count
      repeat = 0 # indictor: 0 => not exist, 1 => exist
      while check <= alibtype - 1:
        if (BAD[libangle[check].atomtypes[0]] == BAD[nlt]) and (BAD[libangle[check].atomtypes[1]] == BAD[nct]) \
and (BAD[libangle[check].atomtypes[2]] == BAD[nrt]) or (((BAD[libangle[check].atomtypes[0]]) == BAD[nrt]) \
and (BAD[libangle[check].atomtypes[1]] == BAD[nct]) and (BAD[libangle[check].atomtypes[2]] == BAD[nlt])):
#        if (libangle[check].atomtypes[0] == nlt) and (libangle[check].atomtypes[1] == nct) \
#and (libangle[check].atomtypes[2] == nrt) or (((libangle[check].atomtypes[0]) == nrt) \
#and (libangle[check].atomtypes[1] == nct) and (libangle[check].atomtypes[2] == nlt)):
          repeat = 1
          ntype = check #angle type that existes and suits to this angle
        check = check + 1  
      # if the angle has been already created in library
      if repeat == 1: 
        nna = Angle()
        nna.atomids = [nl,nc,nr]
        nna.type = ntype
        nnnew[ncount] = nna
        #print ncount, nnnew[ncount].type, nnnew[ncount].atomids[0],nnnew[ncount].atomids[1], \
        #nnnew[ncount].atomids[2]
        ncount = ncount + 1
        continue  
      # if the angle has not yet exist, add it to the library of angle type
      elif repeat == 0:
        nnn = Libangle()
        nnn.atomtypes = [nlt,nct,nrt]
        libangle[alibtype] = nnn
        print alibtype, libangle[alibtype].atomtypes[0], libangle[alibtype].atomtypes[1], \
            libangle[alibtype].atomtypes[2]
        nna = Angle()
        nna.atomids = [nl,nc,nr]
        nna.type = alibtype
        nnnew[ncount] = nna
        print  ncount, nnnew[ncount].type, nnnew[ncount].atomids[0],nnnew[ncount].atomids[1], \
        nnnew[ncount].atomids[2]
        alibtype = alibtype + 1
        ncount = ncount + 1
        repeat = 1
        continue 
    count1 = count1 + 1  

##--------------------------------------------------------------------------'
##------------------------------Create-dihedral--------------------------------'
##--------------------------------------------------------------------------'
#print "____________________________________________________________________\n"
#print "Create new dihedral\n"
b = bnew.copy()
l = liball
ddnew = {} # new library for all dihedral
libdih = {}
ddd = Libdih()
ddd.atomtypes = [0,0,0,0] # dihedral type zero using as the dummy for first checking
libdih[0] = ddd
dlibtype = 1
count0 = 1
dcount =1 # count number of dihedral in data file
for lid in l:
  numl = len(l[lid].atomids)
  count1 = 0
  while count1 < numl : #-1
    d3 = l[lid].atomids[count1]
    d2 = lid
    numl2 = len(l[d3].atomids)
    count2 = 0
    while count2 < numl2:
      if (l[d3].atomids[count2] != d2) and (l[d3].atomids[count2] != d3):
        d4 = l[d3].atomids[count2]
#        count3 = count1 + 1
        count3 = 0
        while count3 < numl:
         d1 = l[lid].atomids[count3]
         if (d1 != d3) and (d1 != d4) and (d1 != d2) and (d2 != d3) and (d2 != d4) and (d3 != d4):
          d1t = a[d1].type
          d2t = a[d2].type
          d3t = a[d3].type
          d4t = a[d4].type
          print "BadID",lid,count0, d1, d2, d3, d4
          #print 'type', d1t, d2t, d3t, d4t
          check = 0 # checking count
          repeat = 0 # indictor: dihedral type, 0 => not exist, 1 => exist
          while check <= dlibtype - 1:
            if ((BAD[libdih[check].atomtypes[0]] == BAD[d1t]) and (BAD[libdih[check].atomtypes[1]] == BAD[d2t]) \
and (BAD[libdih[check].atomtypes[2]] == BAD[d3t]) and (BAD[libdih[check].atomtypes[3]] == BAD[d4t])) or \
((BAD[libdih[check].atomtypes[0]] == BAD[d4t]) and (BAD[libdih[check].atomtypes[1]] == BAD[d3t]) \
and (BAD[libdih[check].atomtypes[2]] == BAD[d2t]) and (BAD[libdih[check].atomtypes[3]] == BAD[d1t])):
              repeat = 1
              dtype = check
            check = check + 1
          # if the dihedral has been already created in library
          if repeat == 1:
            checknew = 1 #index for checking  dihedral set
            repeatnew = 0 #indictor: dihedral set (1-2-3-4) or (4-3-2-1), 0 => not exist, 1 => exist
            # check whether new dihedral is repeated or not
            while checknew <= dcount - 1 : 
              if ((ddnew[checknew].atomids[0] == d1) and (ddnew[checknew].atomids[1] == d2) and\
(ddnew[checknew].atomids[2] == d3) and (ddnew[checknew].atomids[3] == d4)) or \
((ddnew[checknew].atomids[0] == d4) and (ddnew[checknew].atomids[1] == d3) and \
(ddnew[checknew].atomids[2] == d2) and (ddnew[checknew].atomids[3] == d1)):
                repeatnew = 1
              checknew = checknew + 1
            if repeatnew != 1:
              dda = Dihedral()
              dda.atomids = [d1,d2,d3,d4]
              dda.type = dtype
              ddnew[dcount] = dda
              #print dcount, ddnew[dcount].type, ddnew[dcount].atomids[0], ddnew[dcount].atomids[1], \
              #ddnew[dcount].atomids[2], ddnew[dcount].atomids[3]
              dcount = dcount + 1
          # if the dihedral has not yet exist, add it to the library of dihedral type
          elif repeat == 0:
            ddd = Libdih()
            ddd.atomtypes = [d1t,d2t,d3t,d4t]
            libdih[dlibtype] = ddd
            #print dlibtype, libdih[dlibtype].atomtypes[0], libdih[dlibtype].atomtypes[1], \
            #libdih[dlibtype].atomtypes[2], libdih[dlibtype].atomtypes[3]
            dda = Dihedral()
            dda.atomids = [d1,d2,d3,d4]
            dda.type = dlibtype
            ddnew[dcount] = dda
            #print dcount, ddnew[dcount].type, ddnew[dcount].atomids[0], ddnew[dcount].atomids[1], \
            #ddnew[dcount].atomids[2], ddnew[dcount].atomids[3] 
            dlibtype = dlibtype + 1
            dcount = dcount + 1 
            repeat = 1
          count0 = count0 + 1
         count3 = count3 + 1
      count2 = count2 + 1
    count1 = count1 + 1

#Add all bonds, angles, dihedrals for writing output
m.bonds = b
m.angles = nnnew
m.dihedrals = ddnew

#Update all number of each, number of type, and coefficient.
m.nangles = ncount - 1
m.nangletypes =  alibtype - 1
m.ndihedrals = dcount - 1
m.ndihedraltypes = dlibtype -1

#Delete all velocities
m.velocities = {}

#Indicate BAD number for angles and dihedrals in OPLSAA
#Angle coefficients
i = 1 
nat = m.nangletypes
ac = {}
while i <= nat:
  #comment = "#"+BAD[libangle[i].atomtypes[0]]+"-"+BAD[libangle[i].atomtypes[1]]+"-"+BAD[libangle[i].atomtypes[2]]
  #print comment
#  ac[i] = [0,0,comment]
#  print ac[i]
#  i = i + 1
  #check with oplsaa force field parameters list
  j = 1
  while j <= naff:
    if ((BAD[libangle[i].atomtypes[0]] == aff[j][0]) and (BAD[libangle[i].atomtypes[1]] == aff[j][1]) and \
 (BAD[libangle[i].atomtypes[2]] == aff[j][2])) or ((BAD[libangle[i].atomtypes[0]] == aff[j][2]) and \
(BAD[libangle[i].atomtypes[1]] == aff[j][1]) and (BAD[libangle[i].atomtypes[2]] == aff[j][0])):
      fangle = aff[j][3] #angle stiffness
      theta0 = aff[j][4] #theta zero
      ac[i] = [fangle,theta0]
      print i, ac[i]
      i = i + 1
      break
    #if the coefficient cannot be found, it will be assumed to be zero
    if j == naff:
      comment = str(BAD[libangle[i].atomtypes[0]])+"-"+str(BAD[libangle[i].atomtypes[1]])\
+"-"+str(BAD[libangle[i].atomtypes[2]])+' '+str(libangle[i].atomtypes[0])+"-"+str(libangle[i].atomtypes[1])+"-"\
+str(libangle[i].atomtypes[2])
      ac[i] = [0.0,0.0,"#Coefficient-is-not-found",comment]
      print i, ac[i]
      i = i + 1 
    j = j + 1
m.angle_coeffs = ac

#Dihedral coefficients
i = 1
ndt = m.ndihedraltypes
dc = {}
while i <= ndt:
#  comment = "#"+BAD[libdih[i].atomtypes[0]]+"-"+BAD[libdih[i].atomtypes[1]]+"-"+BAD[libdih[i].atomtypes[2]]+"-" \
#+BAD[libdih[i].atomtypes[3]]
#  print comment
#  dc[i] = [0,0,0,0,comment]
#  print dc[i]
#  i = i + 1
  #check with oplsaa force field parameters list
  j = 1
  while j <= ndff:
    if ((BAD[libdih[i].atomtypes[0]] == dff[j][0]) and (BAD[libdih[i].atomtypes[1]] == dff[j][1]) and \
(BAD[libdih[i].atomtypes[2]] == dff[j][2]) and (BAD[libdih[i].atomtypes[3]] == dff[j][3])) or \
((BAD[libdih[i].atomtypes[0]] == dff[j][3]) and (BAD[libdih[i].atomtypes[1]] == dff[j][2]) and \
(BAD[libdih[i].atomtypes[2]] == dff[j][1]) and (BAD[libdih[i].atomtypes[3]] == dff[j][0])):
      r1 = dff[j][4]
      r2 = dff[j][7]
      r3 = dff[j][10]
      r4 = 0
      dc[i] = [r1,r2,r3,r4]
      print i, dc[i]
      i = i + 1
      break
    #if the coefficient cannot be found, it will be assumed to be zero
    if j == ndff:
      dc[i] = [0,0,0,0,"#Coefficient-is-not-found"]
      print i, dc[i]
      i = i + 1
    j = j + 1
m.dihedral_coeffs = dc

#Write the output file
outfile=filename.split(".")[0]
timestr = time.strftime('%d%b%y')
outputfilename = 'complete_'+outfile+'_'+timestr+'.data'
write.moleculefile(outputfilename,m)
