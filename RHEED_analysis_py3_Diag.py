from astropy.io import fits
import tkinter as tk
import tkinter.filedialog as askopenfilename
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage
from scipy import signal
from pylab import figtext
from matplotlib.colors import LogNorm
#import pyperclip


#first we open file
plt.ion()
root=tk.Tk()
fig1=plt.figure(1,figsize=(10,10))
filename = tk.filedialog.askopenfilename(initialdir='/home/will/Documents/School/Data/A383/RHEED_A383/Growth')
hdu_list = fits.open(filename)
image_data = hdu_list[0].data
ax1 = fig1.add_subplot(211)
plt.imshow(image_data, cmap='gray', norm=LogNorm(), aspect='auto')
plt.suptitle(str(filename[len(filename)-7:len(filename)-3]))



coords = []

def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata


    global coords
    coords.append((ix, iy))

    if len(coords) == 2:
        fig1.canvas.mpl_disconnect(cid)
        print(coords)
        
        inc=10000
        d = (np.sqrt((np.abs(coords[0][0]-coords[1][0])**2)+(np.abs(coords[0][1]-coords[1][1])**2)))
        print (d)
        
        
        # to take euclidian change second entry to coords[1][1] not coords[0][1]
        y, x = np.linspace((coords[0][1]), (coords[1][1]), inc), np.linspace((coords[0][0]), (coords[1][0]), inc)

# Zinoise includes background

        zinoise = scipy.ndimage.map_coordinates(image_data, np.vstack((y, x)))

        #zi removes some of the noise

        zi=zinoise[:]-np.min(zinoise)


        line = np.linspace(0,d, num=inc)
        # to plot euclidian line change last entry to coords[1][1] not coords[0][1]
        plt.plot([coords[0][0],coords[1][0]],[coords[0][1],coords[1][1]],'ro-')
        plt.plot([coords[0][0],coords[1][0]],[coords[0][1],coords[0][1]],'ro-')
        plt.plot(coords[0][0], coords[0][1], 'bo')

        plt.plot(coords[1][0], coords[1][1], 'bo')
        
        #Here we find the upper and lower max
        
        search1 = 0
        search2 = 5000
        search3 = 0
        search4 = np.floor(inc/2)
        
        Lmax=np.max(zi[search1:search2])
        Umax=np.max(zi[np.floor(inc//2)+search3:np.floor(inc/2)+search4])

        # Here we find the index of the lower max and upper max

        for Lmaxindex in range(search1, search2):
            if zi[Lmaxindex]==Lmax:
                break

        for Umaxindex in range(int(np.floor(search4+search3)), int(np.floor(search4+search4))):
            if zi[Umaxindex]==Umax:
                break
            
        # here we find the half max of each

        HLmax=np.floor(Lmax/2)
        HUmax=np.floor(Umax/2)

        #Here we find the residuals from the zi with noise removed and the half max for the lower peak and upper peak

        #Lres=zi[0:(inc/2)]-HLmax
        
        Lres=zi[0:inc/2]-HLmax
        Ures=zi[inc/2:inc]-HUmax
        
                 
                 
        #Ures=zi[inc/2:np.floor(inc)]-HUmax
        # Here we find the index of the smallest residual to the left of the lower peak

        # Want to be able to change the search values to specify range, eg. not search middle... doesn't work yet.

        

        for i in range(int(search1),int(Lmaxindex)):
            if np.abs(Lres[i])==np.min(np.abs(Lres[search1:Lmaxindex])):
                break
        #print ('i '+str(i))
        for j in range(int(Lmaxindex), int(search2)):
            if np.abs(Lres[j])==np.min(np.abs(Lres[Lmaxindex:search2])):
                break
        print ('j '+str(j))
        for k in range(int(search3),int(Umaxindex-inc/2)):
            if np.abs(Ures[k])==np.min(np.abs(Ures[search3:Umaxindex-np.floor(inc/2)])):
                print (k)
                break
        #print ('k '+str(k))
        for l in range(int(Umaxindex-np.floor(inc/2)),int(search4)):
            if np.abs(Ures[l])==np.min(np.abs(Ures[Umaxindex-inc/2:search4])):
                break
        print ('l '+str(l))
        
        Lpeakindex=(j+i)/2
        Upeakindex=(l+k)/2+inc/2
        theta=np.arctan((np.abs(coords[1][1]-coords[0][1]))/np.abs(coords[1][0]-coords[0][0]))
        thetadeg=theta*180/3.14159
        
        distance=np.abs(line[Lpeakindex]-line[Upeakindex])
        xdistance=distance*np.cos(theta)
        maxdistance=line[Lmaxindex]-line[Umaxindex]
        ax2 = fig1.add_subplot(212)
        plt.plot(line,zi)
        plt.plot(line[i], zi[i], 'bo')
        plt.plot(line[j], zi[j], 'bo')
        plt.plot(line[int(k+np.floor(inc/2))], zi[int(k+np.floor(inc/2))], 'bo')
        plt.plot(line[l + inc / 2], zi[l + inc / 2], 'bo')
        plt.plot(line[Lpeakindex], zi[Lpeakindex], 'ro')
        plt.plot(line[Upeakindex], zi[Upeakindex], 'ro')
        
        plt.plot(line[Lmaxindex],zi[Lmaxindex],'rx')
        plt.plot(line[Umaxindex],zi[Umaxindex],'rx')
        plt.axvline(line[search3 + inc/2], color='black', linestyle='--')
        plt.axvline(line[search4 + inc / 2 -2], color='black', linestyle='--')

        plt.axvline(line[search1], color='black', linestyle='--')
        plt.axvline(line[search2], color='black', linestyle='--')
        
        figtext(0, 0, '\nLeft peak ' + str(line[Lpeakindex])+"         Left Max "+str(line[Lmaxindex])+ '\nright peak '+str(line[Upeakindex])+"         Right Max"+str(line[Umaxindex])+'\nDistance = '+str(distance)+'         Max Distance '+str(maxdistance)+'            X-Distance ' +str(xdistance))
        print('xdist '+str(xdistance))
        print('theta '+str(thetadeg))
        print(zi[:])
        print(Ures[:])
       # pyperclip.copy(str(distance))
cid = fig1.canvas.mpl_connect('button_press_event', onclick)
plt.show(block=True)
root.destroy()


print('works')


