import os

for i in range(1865,2230):
    filename = "Inside."+'{:04d}'.format(i)
    filename2 = "Inside."+'{:04d}'.format(i-1865)
    #os.system ("mv %s %s" % (filename+'.png', filename+'.tga'))
    os.system ("convert %s %s" % (filename+'.tga', filename2+'.png'))
