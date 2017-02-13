from astropy.io import fits
from tkFileDialog import askopenfilename
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scipy.ndimage
from scipy import signal
from pylab import figtext

#naming figure (1)
fig1=plt.figure(1)




#opening files
filename = askopenfilename()
hdu_list = fits.open(filename)
image_data = hdu_list[0].data
ax1 = fig1.add_subplot(211)
plt.imshow(image_data, cmap='gray', norm=LogNorm())
fig1.suptitle(filename[len(filename)-8:len(filename)-3], fontsize=14, fontweight='bold')

#define varible to place coordinates
coords = []
print (np.size(image_data))

#take clicking coordinates
def onclick(event):
    #global ix, iy
    ix, iy = event.xdata, event.ydata
    print 'x = %d, y = %d '%(
        ix, iy)

    coords.append((ix, iy))
    #print(image_data[(coords[0][0]), (coords[0][1])])
    if len(coords) == 2:

        inc=10000
        d = (np.sqrt((np.abs(coords[0][0]-coords[1][0])**2)+(np.abs(coords[0][1]-coords[1][1])**2)))
        print d
        # to take euclidian change second entry to coords[1][1] not coords[0][1]
        y, x = np.linspace((coords[0][1]), (coords[0][1]), inc), np.linspace((coords[0][0]), (coords[1][0]), inc)

        # Zinoise includes background

        zinoise = scipy.ndimage.map_coordinates(image_data, np.vstack((y, x)))

        #zi removes some of the noise

        zi=zinoise[:]-np.min(zinoise)



        line = np.linspace(0,d, num=inc)
        # to plot euclidian line change last entry to coords[1][1] not coords[0][1]
        plt.plot([coords[0][0],coords[1][0]],[coords[0][1],coords[0][1]],'ro-')
        plt.plot(coords[0][0], coords[0][1], 'bo')

        plt.plot(coords[1][0], coords[1][1], 'bo')

        #Here we find the upper and lower max

        search1 = 0
        search2 = 4000
        search3 = 2500
        search4 = inc/2

        Lmax=np.max(zi[search1:search2])
        Umax=np.max(zi[inc/2+search3:inc/2+search4])

        # Here we find the index of the lower max and upper max

        for Lmaxindex in range(search1, search2):
            if zi[Lmaxindex]==Lmax:
                break

        for Umaxindex in range(inc / 2+search3, inc/2+search4):
            if zi[Umaxindex]==Umax:
                break

        # here we find the half max of each

        HLmax=Lmax/2
        HUmax=Umax/2

        #Here we find the residuals from the zi with noise removed and the half max for the lower peak and upper peak

        Lres=zi[0:inc/2]-HLmax
        Ures=zi[inc/2:inc]-HUmax

        # Here we find the index of the smallest residual to the left of the lower peak

        # we want to zero values in the centre around some point 'null'
        # Want to be able to change the search values to specify range, eg. not search middle... doesn't work yet.



        for i in range(search1,Lmaxindex):
            if Lres[i]==np.min(Lres[search1:Lmaxindex]):
                break

        for j in range(Lmaxindex, search2):
            if Lres[j]==np.min(Lres[Lmaxindex:search2]):
                break

        for k in range(search3,Umaxindex-inc/2):
            if Ures[k]==np.min(Ures[search3:Umaxindex-inc/2]):
                print k
                break

        for l in range(Umaxindex-inc/2,search4):
            if Ures[l]==np.min(Ures[Umaxindex-inc/2:search4]):
                break

        Lpeakindex=(j+i)/2
        Upeakindex=(l+k)/2+inc/2
        distance=line[Lpeakindex]-line[Upeakindex]
        ax2 = fig1.add_subplot(212)
        plt.plot(line,zi)
        #plt.plot(line[inc/2], zi[inc/2], 'bo')
        plt.plot(line[i], zi[i], 'bo')
        plt.plot(line[j], zi[j], 'bo')
        plt.plot(line[k+inc/2], zi[k+inc/2], 'bo')
        plt.plot(line[l + inc / 2], zi[l + inc / 2], 'bo')
        plt.plot(line[Lpeakindex], zi[Lpeakindex], 'ro')
        plt.plot(line[Upeakindex], zi[Upeakindex], 'ro')
        #plt.plot(line[search1], zi[search1],'rx')
        #plt.plot(line[search2], zi[search2], 'rx')
        #plt.plot(line[search3+inc/2], zi[search3+inc/2], 'rx')
        plt.axvline(line[search3 + inc/2], color='k', linestyle='--')
        plt.axvline(line[search4 + inc / 2 -2], color='k', linestyle='--')
        plt.axvline(line[search1], color='k', linestyle='--')
        plt.axvline(line[search2], color='k', linestyle='--')

        figtext(0, 0, '\nLeft peak ' + str(line[Lpeakindex])+ '\nright peak '+str(line[Upeakindex])+'\nDistance = '+str(distance))





cid = fig1.canvas.mpl_connect('button_press_event', onclick)



plt.show()
