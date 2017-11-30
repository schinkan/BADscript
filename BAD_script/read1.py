#!/usr/bin/env python

#Program written by Benjamin Jensen
#Revision 1.0
#January 7th 2013
#Michigan Technological University
#1400 Townsend Dr.
#Houghton, MI 49913

#Edited for taking some comment out from pair_coefficient by Choi, 14 Mar 2017
#Edited for reading OPLSAA forcefield parameter [angle & dihedral] file by Choi, 04 Apr 2017
#Edited for reading OPLSAA forcefield parameter [bond] file by Choi, 05 Apr 2017

class Atom: pass 
class Bond: pass #.type .atomids = [atom1id, atom2id]
class Angle: pass #.type .atomids = [atom1id, atom2id, atom3id]
class Dihedral: pass #.type .atomids = [atom1id, atom2id, atom3id, atom4id]

def strip_comment(line):
  #remove comments
  end = line.find('# ')
  if end >= 0:
      #comment = line[end:]
      line = line[:end]
  return line


class Molecule_File:
  def __init__(self, inmolfile):
    self.atoms = {} #{atom number : atom object}
    self.bonds = {} #{bond number : bond object}
    self.angles = {} #{angle number : angle object}
    self.dihedrals = {}#{dihedral number : dihedral object}

    self.velocities = {} #{atom number : tuple of velocities}

    # Parameters
    self.masses = {} #{atom type : atom mass}
    self.pair_coeffs = {} #{atom type : list of coeffs}
    self.bond_coeffs = {} #{bond type : list of coeffs}  {1: [340,1.5], 2: [450,1.2], ...}
    self.angle_coeffs = {} #{angle type : list of coeffs}
    self.dihedral_coeffs = {} #{dihedral type : list of coeffs}

    self.total_line = ''
    self.ttype_line = ''
    self.xbox_line  = ''
    self.ybox_line  = ''
    self.zbox_line  = ''
    self.extra_lines = ''

    self.total		= 0
    self.natoms		= 0 
    self.natomtypes	= 0
    self.nbonds		= 0
    self.nbondtypes	= 0
    self.nangles	= 0
    self.nangletypes	= 0
    self.ndihedrals	= 0
    self.ndihedraltypes = 0

    self.parsefile(inmolfile)   
    self.inmolfile = inmolfile
 
  def parsefile(self, inmolfile):
    f = open(inmolfile,'r')

    massflag = False
    paircoeff_flag = False
    coeff_flag = False
    atomflag = False
    bondflag = False
    angleflag = False
    dihedralflag = False
    velocityflag = False

    skip = 0
    for line in f:
      #skip comment lines
      skip -= 1
      if skip >= 0: continue
      
      #remove comments
      line = strip_comment(line)
      line = line.strip()
      
      #begining of a section, flag the start and skip one line
      if line == '':
        massflag = False
        paircoeff_flag = False
        coeff_flag = False
        atomflag = False
        bondflag = False
        angleflag = False
        dihedralflag = False
        velocityflag = False
      elif 'atoms' in line:
	self.total_line = line
	self.total = int(line.split()[0])
	self.natoms = self.total
	continue
      elif 'bonds' in line:
	self.nbonds = int(line.split()[0])
	continue
      elif 'angles' in line:
	self.nangles = int(line.split()[0])
	continue
      elif 'dihedrals' in line:
	self.ndihedrals = int(line.split()[0])
	continue
      elif 'atom types' in line:
	self.ttype_line = line
	self.natomtypes = int(line.split()[0])
	continue
      elif 'bond types' in line:
	self.nbondtypes = int(line.split()[0])
	continue
      elif 'angle types' in line:
	self.nangletypes = int(line.split()[0])
	continue
      elif 'dihedral types' in line:
	self.ndihedraltypes = int(line.split()[0])
	continue
      elif 'per atom' in line:
	self.extra_lines += line + '\n'
        continue
      elif 'xlo' in line:
	self.xbox_line = line 
	continue
      elif 'ylo' in line:
	self.ybox_line = line 
	continue
      elif 'zlo' in line:
	self.zbox_line = line 
	continue
      elif line == 'Masses':
        massflag = True
        skip = 1
	continue
      elif line == 'Pair Coeffs':
        paircoeff_flag = True
        coeffs = self.pair_coeffs
        skip = 1
	continue
      elif line == 'Bond Coeffs':
        coeff_flag = True
        coeffs = self.bond_coeffs
        skip = 1
	continue
      elif line == 'Angle Coeffs':
        coeff_flag = True
        coeffs = self.angle_coeffs
        skip = 1
	continue
      elif line == 'Dihedral Coeffs':
        coeff_flag = True
        coeffs = self.dihedral_coeffs
        skip = 1
	continue
      elif line == 'Atoms':
        atomflag = True
        skip = 1
        continue
      elif line == 'Bonds':
        bondflag = True
        skip = 1
        continue
      elif line == 'Angles':
        angleflag = True
        skip = 1
        continue
      elif line == 'Dihedrals':
        dihedralflag = True
        skip = 1
        continue
      elif line == 'Velocities':
        velocityflag = True
        skip = 1
        continue

      if massflag:
        line = line.split()
        id   = int(line[0])
        mass = float(line[1])
        self.masses[id] = mass

      if paircoeff_flag:
        line = line.split()
        id   = int(line[0])
        idcoeffs = []
        for i in line[1:]:
          if '#' not in i:
            idcoeffs.append(float(i))
          else:
            idcoeffs.append(str(i))
        coeffs[id] = idcoeffs

      if coeff_flag:
	line = line.split()
	id   = int(line[0])
	idcoeffs = []
        for i in line[1:]:
	  if '.' not in i:
	    idcoeffs.append(int(i))
	  else:
	    idcoeffs.append(float(i))
        coeffs[id] = idcoeffs

      if atomflag:
        line = line.split()
        id = int(line[0])
        typenum = int(line[1])
        type = int(line[2])
	#charge = float(line[3])
        x = float(line[3])
        y = float(line[4])
        z = float(line[5])
        a = Atom()
        a.typenum = typenum
        a.type = type
	#a.charge = charge
        a.x = x
        a.y = y
        a.z = z
        self.atoms[id] = a
        
      elif bondflag:
        line = line.split()
        id = int(line[0])
        type = int(line[1])
        atom1id = int(line[2])
        atom2id = int(line[3])
        b = Bond()
        b.type = type
        b.atomids = [atom1id, atom2id]
        self.bonds[id] = b
      
      elif angleflag:
        line = line.split()
        id = int(line[0])
        type = int(line[1])
        atom1id = int(line[2])
        atom2id = int(line[3])
        atom3id = int(line[4])
        c = Angle()
        c.type = type
        c.atomids = [atom1id, atom2id, atom3id]
        self.angles[id] = c
        
      elif dihedralflag:
        line = line.split()
        id = int(line[0])
        type = int(line[1])
        atom1id = int(line[2])
        atom2id = int(line[3])
        atom3id = int(line[4])
        atom4id = int(line[5])
        d = Dihedral()
        d.type = type
        d.atomids = [atom1id, atom2id, atom3id, atom4id]
        self.dihedrals[id] = d

      elif velocityflag:
        line = line.split()
        id = int(line[0])
        vx = float(line[1])
        vy = float(line[2])
        vz = float(line[3])
        self.velocities[id] = (vx,vy,vz)

    f.close()
  
class Forcefield:
  def __init__(self, inmolfile):
    self.bondff = {} #{bond parameter list}
    self.angleff = {} #{angle parameter list}
    self.dihedralff = {} #{dihedral parameter list}

    self.nbff = 0 # number of bond parameters
    self.naff = 0 # number of angle parameters
    self.ndff = 0 # number of dihedral parameters

    self.parsefile(inmolfile)   
    self.inmolfile = inmolfile

  def parsefile(self, inmolfile):
    f = open(inmolfile,'r')
    
    bid = 0
    aid = 0
    did = 0
    bcoeffs = {}
    acoeffs = {}
    dcoeffs = {}

    for line in f:
      # bond parameters
      if 'bond ' in line:
        bid =bid + 1
        line = line.split()
        id   = bid
        idcoeffs = []
        for i in line[1:]:
         if '#' not in i:
           if '.' not in i:
             idcoeffs.append(int(i))
           else:
             idcoeffs.append(float(i))
         else:
           pass
        bcoeffs[id] = idcoeffs
      # angle parameters
      if 'angle ' in line:
        aid = aid + 1
        line = line.split()
        id   = aid
        idcoeffs = []
        for i in line[1:]:
         if '#' not in i:
           if '.' not in i:
             idcoeffs.append(int(i))
           else:
             idcoeffs.append(float(i))
         else:
           pass
        acoeffs[id] = idcoeffs
      elif 'torsion ' in line:
        did = did + 1
        line = line.split()
        id   = did
        idcoeffs = []
        for i in line[1:]:
         if '#' not in i: 
           if '.' not in i:
             idcoeffs.append(int(i))
           else:
             idcoeffs.append(float(i))
         else:
           pass
        dcoeffs[id] = idcoeffs
    self.bondff = bcoeffs
    self.angleff = acoeffs
    self.dihedralff = dcoeffs
    self.nbff = bid
    self.naff = aid
    self.ndff = did

    f.close()  
