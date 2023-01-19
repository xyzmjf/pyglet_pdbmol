#!/usr/bin/env python
# Simple code to visualise and rotate molecular structures loaded in PDB format
# now converted python2 to python3
# Simple Pyglet window, but resizable

from pyglet.gl import *
# get keypress definitions etc
from pyglet import window
from pyglet.window import key
from pyglet.window import mouse
from math import sqrt
import random
import sys


# define a function for listing molecule data
def listmol(mol):
	natoms=mol['natoms']
	print('number of atoms = ',natoms)
	#print 'list of coordinates '
	for i in range(1,natoms+1):
		listatom(mol,i)
# end of function

def listatom(mol,i):
	xx=mol[i,'x']
	yy=mol[i,'y']
	zz=mol[i,'z']
	rname=mol[i,'resname']
	rnumber=mol[i,'resnumber']
	aname=mol[i,'atomname']
	print('Atom ',i,aname, rname,rnumber,xx,yy,zz)
# end of function

def distsqatoms(mol,i,j):
	x1=mol1[i,'x']
	y1=mol1[i,'y']
	z1=mol1[i,'z']
	x2=mol1[j,'x']
	y2=mol1[j,'y']
	z2=mol1[j,'z']
	distsq=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
	return(distsq)	
# end of function

def distatoms(mol,i,j):
	distsq=distsqatoms(mol,i,j)
	dist=sqrt(distsq)
	return(dist)
# end of function

# define  bonding elements of mol data structure
# build up data for both i-j and j-i bonds
def makebond(mol,i,j):
	# current numbers of bonds for atoms i and j
	nbondsi=mol[i,'nbonds']
	nbondsj=mol[j,'nbonds']
	nbondedi=1+nbondsi
	nbondedj=1+nbondsj	
	mol[i,'nbonds']=nbondedi
	mol[j,'nbonds']=nbondedj
	mol[i,'bonded',nbondedi]=j
	mol[j,'bonded',nbondedj]=i
# end of function


def drawmol(mol):
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	

# define orthographic projection and window scale in user units
	glMatrixMode(GL_PROJECTION)
	glLoadIdentity()
	glOrtho(-30.0,30.0,-30.0,30.0,-300.0,300.0)

	xcentre=0
	ycentre=0
	zcentre=0
	scale=4
	natoms=mol['natoms']
	#glBegin(GL_LINES)
	for i in range (1,natoms+1):
		x1=mol[i,'x']
		y1=mol[i,'y']
		z1=mol[i,'z']
		xdraw1=float(x1*scale+xcentre)
		ydraw1=float(y1*scale+ycentre)
		zdraw1=float(z1*scale+zcentre)
		#print i,xdraw1,ydraw1,zdraw1
		nbondsi=mol[i,'nbonds']
		#print 
		#listatom(mol,i)
		#print 'Number of bonded atoms =',nbondsi	
		for j in range (1,nbondsi+1):
			atomj=mol[i,'bonded',j]
			#listatom(mol,atomj)
			#print 'distance=',distatoms(mol,i,atomj)
			x2=mol[atomj,'x']
			y2=mol[atomj,'y']
			z2=mol[atomj,'z']
			xdraw2=float(scale*x2+xcentre)
			ydraw2=float(scale*y2+ycentre)
			zdraw2=float(scale*z2+zcentre)
			# white 
			#glColor3f(1.0,1.0,1.0)
			# draw in two stages - half bond colour mode
			xmean=0.5*(xdraw1+xdraw2)
			ymean=0.5*(ydraw1+ydraw2)
			zmean=0.5*(zdraw1+zdraw2)
			red1=mol[i,'red']
			green1=mol[i,'green']
			blue1=mol[i,'blue']
			glColor3f(red1,green1,blue1)
			#pyglet.graphics.draw(2, pyglet.gl.GL_LINES,('v2f', (xdraw1,ydraw1,xdraw2,ydraw2)))
			pyglet.graphics.draw(2, pyglet.gl.GL_LINES,('v3f', (xdraw1,ydraw1,zdraw1,xmean,ymean,zmean)))
			red2=mol[atomj,'red']
			green2=mol[atomj,'green']
			blue2=mol[atomj,'blue']
			glColor3f(red2,green2,blue2)
			pyglet.graphics.draw(2, pyglet.gl.GL_LINES,('v3f', (xmean,ymean,zmean,xdraw2,ydraw2,zdraw2)))
			
	#glEnd()
		
# end of function

def centremol(mol):
	natoms=mol['natoms']
	xsum=0.0
	ysum=0.0
	zsum=0.0
	for i in range (1,natoms+1):
		xx=mol[i,'x']
		yy=mol[i,'y']
		zz=mol[i,'z']
		xsum=xsum+xx
		ysum=ysum+yy
		zsum=zsum+zz
	xmean=xsum/(float(natoms))
	ymean=ysum/(float(natoms))
	zmean=zsum/(float(natoms))
	
	# now subtract mean coordinates
	for i in range (1,natoms+1):
		xx=mol[i,'x']
		yy=mol[i,'y']
		zz=mol[i,'z']
		mol[i,'x']=xx-xmean
		mol[i,'y']=yy-ymean
		mol[i,'z']=zz-zmean
# end of function

def listbonds(mol):
	natoms=mol['natoms']
	for i in range (1,natoms+1):
		listbonds1atom(mol,i)
# end of function

def listbonds1atom(mol,i):
	nbondsi=mol[i,'nbonds']
	print() 
	listatom(mol,i)
	print('Number of bonded atoms =',nbondsi)	
	for j in range (1,nbondsi+1):
		atomj=mol[i,'bonded',j]
		listatom(mol,atomj)
		print('distance=',distatoms(mol,i,atomj))
# end of function



def calcbonds(mol):
# slow calculation of likely bonds
	#cutoff=1.75
	natoms=mol['natoms']
# set initial numbers of bonds to zero
	for i in range(1,natoms+1):
		mol[i,'nbonds']=0
# now calculate bonds based on distance criteria
	for i in range(1,1+natoms):
		for j in range (1,1+natoms):
			if (j>i):			
				distsq=distsqatoms(mol,i,j)
				cutoff=mol[i,'radius']+mol[j,'radius']
				if (distsq<(cutoff*cutoff)):
# now populate bonding elements of data structure
# make i-j and j-i bond
					makebond(mol,i,j)
# end of function


def calc_atom_colours(mol):
	# set red, green, blue colour components for all atoms in molecule
	# store these as elements of mol associative array (dictionary)
	natoms=mol['natoms']
	# set initial colours, will hold for non identified atoms
	for i in range(1,1+natoms):
		mol[i,'red']=1.0
		mol[i,'green']=1.0
		mol[i,'blue']=1.0
	# now set colours for specific atoms
	for i in range(1,1+natoms):
		name=mol[i,'atomname']
		name0=name[0]
		name1=name[1]
		if name0=='H':
			mol[i,'red']=1.0
			mol[i,'green']=1.0
			mol[i,'blue']=1.0
		if name1=='H':
			mol[i,'red']=1.0
			mol[i,'green']=1.0
			mol[i,'blue']=1.0
		if name0=='C':
			mol[i,'red']=0.5
			mol[i,'green']=0.5
			mol[i,'blue']=0.5
		if name0=='O':
			mol[i,'red']=1.0
			mol[i,'green']=0.0
			mol[i,'blue']=0.0
		if name0=='N':
			mol[i,'red']=0.0
			mol[i,'green']=0.0
			mol[i,'blue']=1.0
		if name0=='S':
			mol[i,'red']=1.0
			mol[i,'green']=1.0
			mol[i,'blue']=0.0
		if name0=='P':
			mol[i,'red']=1.0
			mol[i,'green']=0.0
			mol[i,'blue']=1.0
		#print 'first two characters in atom name = ',name0,name1
		#print 'RGB=',mol[i,'red'],mol[i,'green'],mol[i,'blue']


def calc_atom_radius(mol):
	# set atom radii for all atoms in molecule
	natoms=mol['natoms']
	# set initial colours, will hold for non identified atoms
	for i in range(1,1+natoms):
		mol[i,'radius']=0.8
	# now set radii for specific atoms
	for i in range(1,1+natoms):
		name=mol[i,'atomname']
		name0=name[0]
		name1=name[1]
		if name0=='H':
			mol[i,'radius']=0.5
		if name1=='H':
			mol[i,'radius']=0.5
		if name0=='C':
			mol[i,'radius']=0.85
		if name0=='O':
			mol[i,'radius']=0.8
		if name0=='N':
			mol[i,'radius']=0.8
		if name0=='S':
			mol[i,'radius']=1.0
		if name0=='P':
			mol[i,'radius']=1.0
		#print 'first two characters in atom name = ',name0,name1
		#print 'Radius=',mol[i,'radius']

	

# define a function - do not forget the colon 
# parse out atom data from an ATOM or HETATM record in a PDB file
def parse_pdb_line(mol,atomnumber,line):
	#print 'In function parse_pdb_line'
	xstr=line[30:38]
	ystr=line[39:47]
	zstr=line[48:56]
	atmname=line[12:17]
	atmname2=atmname.lstrip()	
	#print atmname,atmname2,xstr, ystr, zstr
	mol[atomnumber,'x']=float(xstr)
	mol[atomnumber,'y']=float(ystr)
	mol[atomnumber,'z']=float(zstr)
	mol[atomnumber,'atomname']=atmname2
	mol[atomnumber,'resname']=line[17:20]
	mol[atomnumber,'resnumber']=line[24:27]
# end of function



# command line arguments
programname=sys.argv[0]
if len(sys.argv)==1:
	print('Usage pdbmol1.py file.pdb')
	sys.exit()

# angle to be used as global variable
xangle=0.0
yangle=0.0
zangle=0.0
	
PDBfilename=sys.argv[1]
print('Program name =',programname)
print('Filename = ',PDBfilename)

# define blank dictionary as molecular data structure
mol1={}

print('If you get this far then should be OK to open file ',PDBfilename)
pdbfile=open(PDBfilename,'r')
pdb_patoms=0
pdb_hetatoms=0

atomcount=0
for line in pdbfile.readlines():
	# print line
	line4=line[:4]
	line6=line[:6]
	# print line4,line6
	if line4=='ATOM':
		atomcount=atomcount+1
		pdb_patoms=pdb_patoms+1
	# call function to atom data from line
		parse_pdb_line(mol1,atomcount,line)
	elif line6=='HETATM': 
		atomcount=atomcount+1
		pdb_hetatoms=pdb_hetatoms+1
	# call function to atom data from line
		parse_pdb_line(mol1,atomcount,line)

mol1['natoms']=pdb_patoms+pdb_hetatoms

print('Number of ATOM records =',pdb_patoms)
print('Number of HETATM records=',pdb_hetatoms)
pdb_natoms=pdb_patoms+pdb_hetatoms

# list molecule data
# print mol1
calc_atom_radius(mol1)
#listmol(mol1)
calcbonds(mol1)


# print molecule bond information
#listbonds(mol1)			

#listatom(mol1,2)
#listatom(mol1,5)
#print distatoms(mol1,2,5)

calc_atom_colours(mol1)



# last code before event loop
#--------------------

# Direct OpenGL commands to this window.
window = pyglet.window.Window(width=500,height=500)

# create graphics context
# see P6 (17 of 101) in programming_guide.pdf
context = window.context
config = context.config
config.double_buffer


# window redraw event
@window.event
def on_draw():
    window.clear()
    drawmol(mol1)	



# window resize event
@window.event
def on_resize(width,height):
	print('The window was resized to %dx%d' % (width, height))	

# window mouse pressed event
@window.event
def on_mouse_press(x, y, button, modifiers):
	print('Mouse button pressed x,y=',x,y)


@window.event
def on_mouse_drag(x,y,dx,dy,buttons,modifiers):
	global xangle
	global yangle
	global zangle	
	if dx < -1:
		#print 'left drag'
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		yangle=float(yangle+2.0)
		glClear(GL_COLOR_BUFFER_BIT)
		glRotatef(xangle,1.0,0.0,0.0)
		glRotatef(yangle,0.0,1.0,0.0)
	elif dx > 1:
		#print 'right drag'
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		yangle=float(yangle-2.0)
		glClear(GL_COLOR_BUFFER_BIT)
		glRotatef(xangle,1.0,0.0,0.0)
		glRotatef(yangle,0.0,1.0,0.0)
	elif dy > 1:
		#print 'Up drag'
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		xangle=float(xangle+2.0)
		glClear(GL_COLOR_BUFFER_BIT)
		glRotatef(xangle,1.0,0.0,0.0)
		glRotatef(yangle,0.0,1.0,0.0)
	elif dy < -1:
		#print 'down drag'
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		xangle=float(xangle-2.0)
		glClear(GL_COLOR_BUFFER_BIT)
		glRotatef(xangle,1.0,0.0,0.0)
		glRotatef(yangle,0.0,1.0,0.0)
		


# key pressed event - if space print something
@window.event
def on_key_press(symbol, modifiers):
	global xangle
	global yangle
	global zangle	
	if symbol == key.ENTER:
		print('ENTER pressed - use this to enter command')    
	elif symbol == key.SPACE:
    		print('The spacebar key was pressed.')



pyglet.app.run()

